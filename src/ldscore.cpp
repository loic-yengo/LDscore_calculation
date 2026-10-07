#include <algorithm>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <pthread.h>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define PACK_DENSITY 4
#define MASK0   3 // 3 << 2 * 0
#define MASK1  12 // 3 << 2 * 1
#define MASK2  48 // 3 << 2 * 2
#define MASK3 192 // 3 << 2 * 3

using namespace std;

namespace {

constexpr int8_t MISSING_GENO = 3;
constexpr double L2_FLOOR = 1e-8;
constexpr int LD_SNP_BLOCK = 32;

void decode_plink(char *output, const char *input, const int lengthInput){
  int i, k;
  char tmp, geno;
  int a1, a2;

  for(i=0;i<lengthInput;++i){
    tmp = input[i];
    k   = PACK_DENSITY * i;
    geno      = (tmp & MASK0);
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno      = (tmp & MASK1) >> 2;
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno      = (tmp & MASK2) >> 4;
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno      = (tmp & MASK3) >> 6;
    a1        = !(geno & 1);
    a2        = !(geno >> 1);
    output[k] = (geno == 1) ? 3 : a1 + a2;
  }
}

string strip_cr(string s){
  if(!s.empty() && s.back()=='\r') s.pop_back();
  return s;
}

vector<string> split_ws(const string &line){
  vector<string> tok;
  stringstream ss(line);
  string t;
  while(ss >> t) tok.push_back(t);
  return tok;
}

string trim_copy(string s){
  size_t a = 0;
  while(a<s.size() && isspace(static_cast<unsigned char>(s[a]))) a++;
  size_t b = s.size();
  while(b>a && isspace(static_cast<unsigned char>(s[b-1]))) b--;
  return s.substr(a, b-a);
}

string to_lower_copy(string s){
  for(char &c : s) c = static_cast<char>(tolower(static_cast<unsigned char>(c)));
  return s;
}

string norm_chr(string s){
  s = trim_copy(s);
  string low = to_lower_copy(s);
  if(low.size()>=3 && low.compare(0,3,"chr")==0) low = low.substr(3);
  if(low.empty()) return "";
  bool digits = true;
  for(char c : low){
    if(!isdigit(static_cast<unsigned char>(c))){ digits = false; break; }
  }
  if(digits) return to_string(atoi(low.c_str()));
  if(low=="x") return "X";
  if(low=="y") return "Y";
  if(low=="m" || low=="mt" || low=="mito") return "MT";
  return low;
}

bool is_autosome_key(const string &key){
  if(key.empty()) return false;
  for(char c : key){
    if(!isdigit(static_cast<unsigned char>(c))) return false;
  }
  const int n = atoi(key.c_str());
  return n>=1 && n<=22 && to_string(n)==key;
}

bool parse_full_double(const string &s, double &out){
  if(s.empty()) return false;
  char *end = nullptr;
  out = strtod(s.c_str(), &end);
  return end != s.c_str() && *end == '\0' && isfinite(out);
}

bool parse_full_u32(const string &s, uint32_t &out){
  if(s.empty() || s[0]=='-' || s[0]=='+') return false;
  errno = 0;
  char *end = nullptr;
  const unsigned long v = strtoul(s.c_str(), &end, 10);
  if(end == s.c_str() || *end != '\0' || errno != 0) return false;
  if(v > static_cast<unsigned long>(numeric_limits<uint32_t>::max())) return false;
  out = static_cast<uint32_t>(v);
  return true;
}

bool parse_full_bp(const string &s, int &out){
  if(s.empty()) return false;
  errno = 0;
  char *end = nullptr;
  const long v = strtol(s.c_str(), &end, 10);
  if(end == s.c_str() || *end != '\0' || errno != 0) return false;
  if(v < 0 || v > static_cast<long>(numeric_limits<int>::max())) return false;
  out = static_cast<int>(v);
  return true;
}

string mem_line(size_t n_people, size_t n_snps, size_t bytes_each){
  const double bytes = static_cast<double>(n_people) * static_cast<double>(n_snps) * static_cast<double>(bytes_each);
  ostringstream oss;
  oss.setf(ios::fixed);
  oss<<setprecision(2)<<(bytes / 1e9);
  return to_string(n_people)+" individuals x "+to_string(n_snps)+" SNPs x "+to_string(bytes_each)+
         " byte"+((bytes_each==1)?"":"s")+" = "+oss.str()+" GB";
}

struct SayFn {
  ofstream *log = nullptr;
  bool verbose = true;
  void operator()(const string &msg) const {
    if(log) (*log)<<msg<<endl;
    if(verbose) cout<<msg<<endl;
  }
};

struct Snp {
  string chr;
  string id;
  int bp = 0;
};

struct MapPoint {
  int bp = 0;
  double cm = 0.0;
};

struct BimSet {
  string prefix;
  string bedfile;
  vector<Snp> snps;
  unordered_map<string, vector<int>> chr_snps;
  unordered_map<string, string> chr_label;
  vector<string> dropped;
};

string join_keys(const vector<string> &v){
  string s;
  for(size_t i=0;i<v.size();++i){
    if(i) s += ", ";
    s += v[i];
  }
  return s;
}

double interpolate_cm(int bp, const vector<MapPoint> &map){
  if(map.empty()) return numeric_limits<double>::quiet_NaN();
  if(bp <= map.front().bp) return map.front().cm;
  if(bp >= map.back().bp) return map.back().cm;
  auto it = upper_bound(map.begin(), map.end(), bp,
                        [](int x, const MapPoint &p){ return x < p.bp; });
  const MapPoint &r = *it;
  const MapPoint &l = *(it - 1);
  if(r.bp == l.bp) return l.cm;
  double t = (static_cast<double>(bp) - l.bp) / (static_cast<double>(r.bp) - l.bp);
  return l.cm + t * (r.cm - l.cm);
}

// Boundary points are (Begin, cM at Begin), then each (End, cM at End).
// With interpolation, bp outside the span is clamped to the nearest end.
// Without it, bp in [Begin, End) takes the cM at Begin, the last End takes
// the last cM, and bp outside the span is rejected.
bool lookup_cm(int bp, const vector<MapPoint> &map, bool interpolate, double &cm){
  if(map.size() < 2) return false;
  if(interpolate){
    cm = interpolate_cm(bp, map);
    return isfinite(cm);
  }
  if(bp < map.front().bp || bp > map.back().bp) return false;
  auto it = lower_bound(map.begin(), map.end(), bp,
                        [](const MapPoint &p, int x){ return p.bp < x; });
  if(it == map.end()) return false;
  if(it->bp == bp){
    cm = it->cm;
    return true;
  }
  if(it == map.begin()) return false;
  cm = (it - 1)->cm;
  return true;
}

int map_status(int bp, const vector<MapPoint> &map){
  // 0 = exact, 1 = interpolated, 2 = clamped
  if(map.empty()) return 2;
  if(bp < map.front().bp || bp > map.back().bp) return 2;
  auto it = lower_bound(map.begin(), map.end(), bp,
                        [](const MapPoint &p, int x){ return p.bp < x; });
  if(it != map.end() && it->bp == bp) return 0;
  return 1;
}

struct RawInterval {
  int begin = 0;
  int end = 0;
  double rate = 0.0;
  double cm = 0.0;
};

bool is_int_token(const string &s){
  if(s.empty()) return false;
  size_t i = 0;
  if(s[0]=='+') i = 1;
  if(i>=s.size()) return false;
  for(; i<s.size(); ++i){
    if(!isdigit(static_cast<unsigned char>(s[i]))) return false;
  }
  return true;
}

bool is_float_token(const string &s){
  if(s.empty()) return false;
  char *end = nullptr;
  strtod(s.c_str(), &end);
  return end != s.c_str() && end && *end=='\0';
}

// Five columns, in order: chromosome, begin bp, end bp, cM/Mb, cumulative cM at end.
bool is_interval_row(const vector<string> &tok){
  return tok.size()==5 && is_int_token(tok[1]) && is_int_token(tok[2]) &&
         is_float_token(tok[3]) && is_float_token(tok[4]);
}

// One file holds every autosome. Column 5 is the cumulative cM at End.
// Boundary points stored per chromosome: (first Begin, 0) then each (End, cM).
// The first non-empty line is a header when it is not itself a data row.
bool read_interval_map(const string &mapFile,
                       unordered_map<string, vector<MapPoint>> &by_chr,
                       vector<string> &skipped,
                       string &err){
  by_chr.clear();
  skipped.clear();
  ifstream in(mapFile.c_str());
  if(!in){
    err = "Error reading "+mapFile;
    return false;
  }
  unordered_map<string, vector<RawInterval>> raw;
  unordered_set<string> closed;
  unordered_set<string> skip_seen;
  string prev_key;
  string line;
  int lineno = 0;
  bool header_skipped = false;
  int header_width = 0;
  while(getline(in, line)){
    ++lineno;
    line = trim_copy(strip_cr(line));
    if(line.empty()) continue;
    const vector<string> tok = split_ws(line);
    if(!header_skipped && raw.empty() && !is_interval_row(tok)){
      header_skipped = true;
      header_width = static_cast<int>(tok.size());
      continue;
    }
    if(tok.size()!=5){
      err = mapFile+" line "+to_string(lineno)+" must have 5 columns: chromosome, begin, end, cM/Mb, cM.";
      return false;
    }
    if(!is_interval_row(tok)){
      err = mapFile+" line "+to_string(lineno)+" must have integer begin and end, and numeric cM/Mb and cM.";
      return false;
    }
    const string key = norm_chr(tok[0]);
    if(key.empty()){
      err = mapFile+" line "+to_string(lineno)+": could not parse chromosome \""+tok[0]+"\".";
      return false;
    }
    if(!is_autosome_key(key)){
      if(!skip_seen.count(key)){
        skip_seen.insert(key);
        skipped.push_back(key);
      }
      prev_key.clear();
      continue;
    }
    if(closed.count(key)){
      err = "Chromosome "+key+" appears in more than one block in "+mapFile+".";
      return false;
    }
    if(!prev_key.empty() && prev_key!=key) closed.insert(prev_key);
    RawInterval iv;
    iv.begin = atoi(tok[1].c_str());
    iv.end = atoi(tok[2].c_str());
    iv.rate = atof(tok[3].c_str());
    iv.cm = atof(tok[4].c_str());
    if(iv.begin<0 || iv.end<=iv.begin){
      err = mapFile+" line "+to_string(lineno)+": Begin must be >= 0 and End must be greater than Begin.";
      return false;
    }
    if(!isfinite(iv.rate) || iv.rate<0.0){
      err = mapFile+" line "+to_string(lineno)+": cM/Mb cannot be missing or negative.";
      return false;
    }
    if(isnan(iv.cm) || iv.cm<0.0){
      err = mapFile+" line "+to_string(lineno)+": cM cannot be missing or negative.";
      return false;
    }
    raw[key].push_back(iv);
    prev_key = key;
  }
  if(raw.empty()){
    if(header_skipped && header_width!=5){
      err = mapFile+" must have 5 columns: chromosome, begin, end, cM/Mb, cM.";
    }else{
      err = "No autosome intervals (chromosomes 1-22) in "+mapFile;
    }
    return false;
  }
  for(auto &kv : raw){
    const vector<RawInterval> &ivs = kv.second;
    vector<MapPoint> pts;
    pts.reserve(ivs.size()+1);
    MapPoint left;
    left.bp = ivs.front().begin;
    left.cm = 0.0;
    pts.push_back(left);
    for(size_t i=0;i<ivs.size();++i){
      if(i>0 && ivs[i].begin != ivs[i-1].end){
        err = "Chromosome "+kv.first+" in "+mapFile+" has a gap or overlap between intervals.";
        return false;
      }
      if(i>0 && ivs[i].cm < ivs[i-1].cm){
        err = "Chromosome "+kv.first+" in "+mapFile+" has a decreasing cM column.";
        return false;
      }
      MapPoint right;
      right.bp = ivs[i].end;
      right.cm = ivs[i].cm;
      pts.push_back(right);
    }
    if(pts.back().cm > 5000.0){
      err = "Chromosome "+kv.first+" in "+mapFile+" has a maximum cM above 5000.";
      return false;
    }
    by_chr[kv.first] = std::move(pts);
  }
  return true;
}

bool load_fam(const string &famfile, vector<string> &fid, vector<string> &iid, int &n_dup, string &err){
  ifstream in(famfile.c_str());
  if(!in){
    err = "Error reading "+famfile;
    return false;
  }
  string line;
  int n = -1;
  while(in){
    getline(in, line);
    n++;
  }
  in.close();
  fid.assign(n, "");
  iid.assign(n, "");
  n_dup = 0;
  unordered_map<string,int> seen;
  in.open(famfile.c_str());
  if(!in){
    err = "Error reading "+famfile;
    return false;
  }
  for(int i=0;i<n;++i){
    if(!getline(in, line)){
      err = "Error reading "+famfile;
      return false;
    }
    line = strip_cr(line);
    stringstream ss(line);
    ss >> fid[i] >> iid[i];
    if(ss.fail() || fid[i].empty() || iid[i].empty()){
      err = famfile+" line "+to_string(i+1)+" must start with FID and IID.";
      return false;
    }
    if(seen.find(iid[i])!=seen.end()) n_dup++;
    else seen[iid[i]] = i;
  }
  return true;
}

bool read_mbfile(const string &path, vector<string> &prefixes, string &err){
  ifstream in(path.c_str());
  if(!in){
    err = "Error reading "+path;
    return false;
  }
  string line;
  int lineno = 0;
  unordered_set<string> seen;
  while(getline(in, line)){
    ++lineno;
    line = trim_copy(strip_cr(line));
    if(line.empty() || line[0]=='#') continue;
    if(!seen.insert(line).second){
      err = path+" line "+to_string(lineno)+": duplicate prefix \""+line+"\".";
      return false;
    }
    prefixes.push_back(line);
  }
  if(prefixes.empty()){
    err = "No PLINK prefixes in "+path+".";
    return false;
  }
  return true;
}

// FAM row indices used for allele frequency, missingness, and LD scores.
// An empty keepFile uses every row. Otherwise the file lists FID and IID.
bool select_samples(const vector<string> &fid, const vector<string> &iid,
                    const string &keepFile, vector<int> &sample_ids,
                    vector<string> &notes, string &err){
  sample_ids.clear();
  notes.clear();
  const int n = static_cast<int>(iid.size());
  if(n<=0 || static_cast<int>(fid.size())!=n){
    err = "FAM contains no individuals.";
    return false;
  }
  if(keepFile.empty()){
    sample_ids.resize(n);
    iota(sample_ids.begin(), sample_ids.end(), 0);
    return true;
  }

  unordered_map<string,int> row;
  unordered_map<string,int> copies;
  row.reserve(static_cast<size_t>(n));
  for(int i=0;i<n;++i){
    const string key = fid[i]+"\t"+iid[i];
    copies[key]++;
    if(!row.count(key)) row[key] = i;
  }

  ifstream in(keepFile.c_str());
  if(!in){
    err = "Error reading "+keepFile;
    return false;
  }
  struct IdLine { int lineno; string fid; string iid; };
  vector<IdLine> rows;
  string line;
  int lineno = 0;
  while(getline(in, line)){
    ++lineno;
    line = trim_copy(strip_cr(line));
    if(line.empty() || line[0]=='#') continue;
    const vector<string> tok = split_ws(line);
    if(tok.size()!=2){
      err = keepFile+" line "+to_string(lineno)+" must have 2 columns: FID and IID.";
      return false;
    }
    rows.push_back(IdLine{lineno, tok[0], tok[1]});
  }
  if(rows.empty()){
    err = "No individuals in "+keepFile+".";
    return false;
  }

  auto fid_iid_header = [](const string &a0, const string &b0){
    string a = to_lower_copy(a0);
    string b = to_lower_copy(b0);
    if(!a.empty() && a[0]=='#') a = a.substr(1);
    return a=="fid" && b=="iid";
  };
  int start = 0;
  const string key0 = rows[0].fid+"\t"+rows[0].iid;
  if(fid_iid_header(rows[0].fid, rows[0].iid)){
    start = 1;
    notes.push_back("# [Note] "+keepFile+" line "+to_string(rows[0].lineno)+
                    " treated as a header ("+rows[0].fid+" "+rows[0].iid+").");
  }else if(!row.count(key0)){
    if(rows.size()==1){
      err = keepFile+" line "+to_string(rows[0].lineno)+": "+rows[0].fid+" "+rows[0].iid+
            " is not in the FAM.";
      return false;
    }
    start = 1;
    notes.push_back("# [Note] "+keepFile+" line "+to_string(rows[0].lineno)+
                    " treated as a header ("+rows[0].fid+" "+rows[0].iid+").");
  }
  if(start>=static_cast<int>(rows.size())){
    err = "No individuals in "+keepFile+".";
    return false;
  }

  unordered_set<string> seen_keep;
  int n_dup_rows = 0;
  for(int i=start;i<static_cast<int>(rows.size());++i){
    const string key = rows[i].fid+"\t"+rows[i].iid;
    if(!seen_keep.insert(key).second){
      err = keepFile+" line "+to_string(rows[i].lineno)+": duplicate individual "+
            rows[i].fid+" "+rows[i].iid+".";
      return false;
    }
    const auto it = row.find(key);
    if(it==row.end()){
      err = keepFile+" line "+to_string(rows[i].lineno)+": "+rows[i].fid+" "+rows[i].iid+
            " is not in the FAM.";
      return false;
    }
    if(copies[key]>1) n_dup_rows++;
    sample_ids.push_back(it->second);
  }
  sort(sample_ids.begin(), sample_ids.end());
  if(n_dup_rows>0){
    notes.push_back("# [Warning] "+to_string(n_dup_rows)+
                    " --keep individual"+(n_dup_rows==1 ? "" : "s")+
                    " match more than one FAM row; the first row is used.");
  }
  return true;
}

// Draw nsample eligible FAM rows without replacement. The same draw is used
// for every chromosome. When the pool is no larger than nsample, everyone is kept.
vector<int> draw_ld_sample(const vector<int> &eligible, uint32_t nsample, uint32_t seed, bool &drew){
  drew = false;
  if(eligible.size() <= static_cast<size_t>(nsample)) return eligible;
  vector<int> order = eligible;
  mt19937 rng(seed);
  shuffle(order.begin(), order.end(), rng);
  order.resize(nsample);
  sort(order.begin(), order.end());
  drew = true;
  return order;
}

bool load_bim(const string &bimfile, vector<Snp> &snps, string &err){
  ifstream in(bimfile.c_str());
  if(!in){
    err = "Error reading "+bimfile;
    return false;
  }
  string line;
  int p = -1;
  while(in){
    getline(in, line);
    p++;
  }
  in.close();
  snps.assign(p, Snp());
  in.open(bimfile.c_str());
  if(!in){
    err = "Error reading "+bimfile;
    return false;
  }
  for(int j=0;j<p;++j){
    if(!getline(in, line)){
      err = "Error reading "+bimfile;
      return false;
    }
    line = strip_cr(line);
    stringstream ss(line);
    string cm_ignored, bp_s, a1, a2;
    ss >> snps[j].chr >> snps[j].id >> cm_ignored >> bp_s >> a1 >> a2;
    if(ss.fail() || snps[j].chr.empty() || snps[j].id.empty() || !parse_full_bp(bp_s, snps[j].bp)){
      err = "Failed to parse 6-column BIM line: "+trim_copy(line);
      return false;
    }
  }
  return true;
}

bool index_bim_autosomes(const string &bimfile, const vector<Snp> &snps, BimSet &bim, string &err){
  unordered_set<string> dropped_seen;
  for(int j=0;j<static_cast<int>(snps.size());++j){
    const string key = norm_chr(snps[j].chr);
    if(key.empty()){
      err = "Could not parse chromosome from BIM SNP "+snps[j].id+" ("+snps[j].chr+") in "+bimfile+".";
      return false;
    }
    if(!is_autosome_key(key)){
      if(!dropped_seen.count(key)){
        dropped_seen.insert(key);
        bim.dropped.push_back(key);
      }
      continue;
    }
    if(bim.chr_snps[key].empty()){
      bim.chr_label[key] = snps[j].chr;
    }else{
      const int prev = bim.chr_snps[key].back();
      if(snps[j].bp < snps[prev].bp){
        err = "BIM base-pair positions are not non-decreasing on chromosome "+snps[j].chr+" in "+bimfile+".";
        return false;
      }
    }
    bim.chr_snps[key].push_back(j);
  }
  return true;
}

bool open_bed_file(ifstream &in, const string &bedfile, string &err){
  in.open(bedfile.c_str(), ios::in | ios::binary);
  if(!in){
    err = "Error reading "+bedfile;
    return false;
  }
  unsigned char magic[3] = {0, 0, 0};
  in.read(reinterpret_cast<char*>(magic), 3);
  if(!in || magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01){
    err = bedfile+" is not a SNP-major PLINK BED file.";
    return false;
  }
  return true;
}

bool check_bed_size(const string &bedfile, size_t n_snps, int numBytes, string &err){
  ifstream in;
  if(!open_bed_file(in, bedfile, err)) return false;
  in.seekg(0, ios::end);
  if(!in){
    err = "Error reading "+bedfile;
    return false;
  }
  const streamoff sz = in.tellg();
  const streamoff expect = static_cast<streamoff>(3) +
                           static_cast<streamoff>(n_snps) * static_cast<streamoff>(numBytes);
  if(sz != expect){
    err = bedfile+" has "+to_string(static_cast<long long>(sz))+
          " bytes; expected "+to_string(static_cast<long long>(expect))+
          " for "+to_string(n_snps)+" SNPs.";
    return false;
  }
  return true;
}

bool read_bed_snp(ifstream &in, char *packed, char *unpacked, int numBytes, int j, int &last_j){
  if(last_j<0 || j!=last_j+1){
    const streamoff off = static_cast<streamoff>(3) +
                          static_cast<streamoff>(j) * static_cast<streamoff>(numBytes);
    in.seekg(off, ios::beg);
    if(!in) return false;
  }
  in.read(packed, sizeof(char)*numBytes);
  if(!in) return false;
  decode_plink(unpacked, packed, numBytes);
  last_j = j;
  return true;
}

// Unbiased r^2 from sibIBD: r^2 - (1-r^2)/(n-2). Negative estimates are kept.
double r2_est(const int8_t *g1, const int8_t *g2, int n){
  double s1=0, s2=0, nobs=0;
  for(int i=0;i<n;++i){
    if(g1[i]==MISSING_GENO || g2[i]==MISSING_GENO) continue;
    s1 += g1[i];
    s2 += g2[i];
    nobs += 1.0;
  }
  if(nobs < 3.0) return 0.0;
  const double m1 = s1 / nobs;
  const double m2 = s2 / nobs;
  double v1=0, v2=0, cv=0;
  for(int i=0;i<n;++i){
    if(g1[i]==MISSING_GENO || g2[i]==MISSING_GENO) continue;
    const double d1 = g1[i] - m1;
    const double d2 = g2[i] - m2;
    v1 += d1*d1;
    v2 += d2*d2;
    cv += d1*d2;
  }
  if(v1 <= 0.0 || v2 <= 0.0) return 0.0;
  const double r = cv / sqrt(v1 * v2);
  const double r2 = r * r;
  return r2 - (1.0 - r2) / (nobs - 2.0);
}

// Pair contributions for target SNPs [begin, end). Each SNP's own 1 is added by the caller.
void add_window_pairs(vector<double> &ld, const int8_t *gts, int n,
                      const double *cm, double window_cm, int begin, int end){
  for(int i=begin; i<end; ++i){
    const int8_t *gi = gts + static_cast<size_t>(i)*n;
    for(int j=i-1; j>=0; --j){
      const double dist = cm[i] - cm[j];
      if(dist >= window_cm) break;
      const double r2 = r2_est(gi, gts + static_cast<size_t>(j)*n, n);
      ld[i] += r2;
      ld[j] += r2;
    }
  }
}

// Each SNP starts at 1. A pair contributes once cm[i] - cm[j] is below window_cm.
vector<double> compute_ld_scores(const int8_t *gts, int n, const vector<double> &cm, double window_cm){
  const int M = static_cast<int>(cm.size());
  vector<double> ld(M, 1.0);
  if(M>0) add_window_pairs(ld, gts, n, cm.data(), window_cm, 0, M);
  return ld;
}

struct LdWork {
  const int8_t *gts = nullptr;
  const double *cm = nullptr;
  vector<double> *partial = nullptr;
  double window_cm = 0.0;
  int n = 0;
  int M = 0;
  int nblocks = 0;
  int nthr = 0;
};

struct LdWorkArg {
  LdWork *work = nullptr;
  int tid = 0;
};

void *ld_worker(void *arg){
  auto *wa = static_cast<LdWorkArg*>(arg);
  LdWork &w = *wa->work;
  vector<double> &ld = w.partial[wa->tid];
  for(int b=wa->tid; b<w.nblocks; b+=w.nthr){
    const int begin = b * LD_SNP_BLOCK;
    int end = begin + LD_SNP_BLOCK;
    if(end > w.M) end = w.M;
    add_window_pairs(ld, w.gts, w.n, w.cm, w.window_cm, begin, end);
  }
  return nullptr;
}

// Workers share the genotype matrix and each writes into its own score vector.
// Block b is owned by thread b % nthread, so the same thread count repeats.
bool compute_ld_scores_threaded(const int8_t *gts, int n, const vector<double> &cm,
                                double window_cm, int nthread,
                                vector<double> &ld, int &n_used, string &err){
  const int M = static_cast<int>(cm.size());
  n_used = 1;
  if(nthread<=1 || M<=1){
    ld = compute_ld_scores(gts, n, cm, window_cm);
    return true;
  }
  const int nblocks = (M + LD_SNP_BLOCK - 1) / LD_SNP_BLOCK;
  int nthr = nthread;
  if(nthr > nblocks) nthr = nblocks;
  if(nthr<=1){
    ld = compute_ld_scores(gts, n, cm, window_cm);
    return true;
  }
  n_used = nthr;
  vector<vector<double>> partial(nthr, vector<double>(M, 0.0));
  LdWork work;
  work.gts = gts;
  work.cm = cm.data();
  work.partial = partial.data();
  work.window_cm = window_cm;
  work.n = n;
  work.M = M;
  work.nblocks = nblocks;
  work.nthr = nthr;
  vector<LdWorkArg> args(nthr);
  vector<pthread_t> threads(nthr);
  int started = 0;
  bool ok = true;
  for(int t=0; t<nthr; ++t){
    args[t].work = &work;
    args[t].tid = t;
    if(pthread_create(&threads[t], nullptr, ld_worker, &args[t])!=0){
      ok = false;
      err = "pthread_create failed.";
      break;
    }
    started++;
  }
  for(int t=0; t<started; ++t) pthread_join(threads[t], nullptr);
  if(!ok) return false;
  ld.assign(M, 1.0);
  for(int t=0; t<nthr; ++t){
    const vector<double> &part = partial[t];
    for(int i=0; i<M; ++i) ld[i] += part[i];
  }
  return true;
}

struct ChrResult {
  vector<Snp> snps;
  vector<double> cm;
  vector<int8_t> gts;
  vector<double> l2;
  int n_maf = 0;
  int n_missed = 0;
  int n_nomap = 0;
  int n_exact = 0;
  int n_interp = 0;
  int n_clamp = 0;
  int n_floor = 0;
};

bool score_chromosome(const string &bedfile,
                      const vector<Snp> &snps,
                      const vector<int> &snp_idx,
                      const vector<MapPoint> *gmap,
                      const vector<int> &filter_ids,
                      const vector<int> &ld_ids,
                      bool do_interpolate,
                      double min_maf, double max_missing,
                      int numBytes, char *packed, char *unpacked,
                      ChrResult &out, string &err){
  out = ChrResult();
  ifstream influx;
  if(!open_bed_file(influx, bedfile, err)) return false;
  int last_j = -1;
  const int n_filter = static_cast<int>(filter_ids.size());
  const int n_ld = static_cast<int>(ld_ids.size());
  out.snps.reserve(snp_idx.size());
  out.cm.reserve(snp_idx.size());
  out.gts.reserve(snp_idx.size() * static_cast<size_t>(max(n_ld, 1)));
  unordered_set<string> seen_id;
  auto tally = [](int st, int &n_exact, int &n_interp, int &n_clamp){
    if(st==0) n_exact++;
    else if(st==1) n_interp++;
    else n_clamp++;
  };
  for(int j : snp_idx){
    if(!read_bed_snp(influx, packed, unpacked, numBytes, j, last_j)){
      err = "Failed reading SNP "+to_string(j)+" from "+bedfile;
      return false;
    }
    int nmiss = 0, nobs = 0;
    double sx = 0.0;
    for(int s=0;s<n_filter;++s){
      const int x = static_cast<unsigned char>(unpacked[filter_ids[s]]);
      if(x==MISSING_GENO) nmiss++;
      else{ nobs++; sx += x; }
    }
    const int ntot = nobs + nmiss;
    const double miss_frac = ntot>0 ? static_cast<double>(nmiss) / ntot : 1.0;
    if(miss_frac > max_missing || nobs==0){ out.n_missed++; continue; }
    const double p = (sx / nobs) / 2.0;
    const double maf = min(p, 1.0-p);
    if(maf < min_maf){ out.n_maf++; continue; }
    double pos = 0.0;
    if(gmap){
      if(!lookup_cm(snps[j].bp, *gmap, do_interpolate, pos)){ out.n_nomap++; continue; }
    }else{
      pos = static_cast<double>(snps[j].bp);
    }
    if(!seen_id.insert(snps[j].id).second){
      err = "Duplicate SNP id \""+snps[j].id+"\" on chromosome "+snps[j].chr+".";
      return false;
    }
    if(gmap) tally(map_status(snps[j].bp, *gmap), out.n_exact, out.n_interp, out.n_clamp);
    out.snps.push_back(snps[j]);
    out.cm.push_back(pos);
    for(int s=0;s<n_ld;++s){
      out.gts.push_back(static_cast<int8_t>(static_cast<unsigned char>(unpacked[ld_ids[s]])));
    }
  }
  for(int k=1;k<static_cast<int>(out.snps.size());++k){
    if(out.cm[k] < out.cm[k-1]){
      err = gmap
        ? "Genetic map position is not non-decreasing on chromosome "+out.snps[k].chr+". Check the map and BIM order."
        : "Base-pair position is not non-decreasing on chromosome "+out.snps[k].chr+". Check the BIM order.";
      return false;
    }
  }
  return true;
}

bool read_snp_list(const string &path, unordered_set<string> &ids, int &n_lines, string &err){
  ifstream in(path.c_str());
  if(!in){
    err = "Error reading "+path;
    return false;
  }
  string line;
  while(getline(in, line)){
    line = trim_copy(strip_cr(line));
    if(line.empty() || line[0]=='#') continue;
    ids.insert(line);
    n_lines++;
  }
  return true;
}

bool read_extract_list(const string &path, vector<string> &files, string &err){
  ifstream in(path.c_str());
  if(!in){
    err = "Error reading "+path;
    return false;
  }
  string line;
  int lineno = 0;
  unordered_set<string> seen;
  while(getline(in, line)){
    ++lineno;
    line = trim_copy(strip_cr(line));
    if(line.empty() || line[0]=='#') continue;
    if(!seen.insert(line).second){
      err = path+" line "+to_string(lineno)+": duplicate SNP-list path \""+line+"\".";
      return false;
    }
    ifstream probe(line.c_str());
    if(!probe){
      err = "Error reading "+line;
      return false;
    }
    files.push_back(line);
  }
  if(files.empty()){
    err = "No SNP-list files in "+path+".";
    return false;
  }
  return true;
}

struct ExtractCounts {
  int n_kept = 0;
  int n_dropped = 0;
  int n_absent = 0;
};

void apply_extract(vector<BimSet> &bims, const unordered_set<string> &ids,
                   ExtractCounts &counts, unordered_set<string> &emptied){
  unordered_set<string> bim_ids;
  unordered_set<string> had, still;
  counts = ExtractCounts();
  emptied.clear();
  for(auto &bim : bims){
    for(const auto &snp : bim.snps) bim_ids.insert(snp.id);
    vector<string> drop_keys;
    for(auto &kv : bim.chr_snps){
      had.insert(kv.first);
      vector<int> keep;
      keep.reserve(kv.second.size());
      for(int j : kv.second){
        if(ids.count(bim.snps[j].id)){
          keep.push_back(j);
          counts.n_kept++;
        }else{
          counts.n_dropped++;
        }
      }
      if(keep.empty()){
        drop_keys.push_back(kv.first);
      }else{
        kv.second.swap(keep);
        still.insert(kv.first);
      }
    }
    for(const string &key : drop_keys) bim.chr_snps.erase(key);
  }
  for(const string &id : ids){
    if(!bim_ids.count(id)) counts.n_absent++;
  }
  for(const string &key : had){
    if(!still.count(key)) emptied.insert(key);
  }
}

void print_help(){
  cerr<<"\t--bfile        : One binary PLINK prefix ([prefix].bed/.bim/.fam). A merged file or one chromosome."<<endl;
  cerr<<"\t--mbfile       : Text file of PLINK prefixes, one per line. Each prefix is a chromosome or a chunk."<<endl;
  cerr<<"\t                 Blank lines and lines starting with # are ignored. Use either --bfile or --mbfile."<<endl;
  cerr<<"\t                 Every FAM must list the same individuals in the same order."<<endl;
  cerr<<"\t                 A chromosome split across prefixes is scored separately in each prefix."<<endl;
  cerr<<"\t--keep         : Individuals used for allele frequency and missingness."<<endl;
  cerr<<"\t                 Two columns, FID and IID, separated by space or tab."<<endl;
  cerr<<"\t                 Blank lines and lines starting with # are ignored."<<endl;
  cerr<<"\t                 A header is optional. The first line is a header when it is FID and IID,"<<endl;
  cerr<<"\t                 or when it is not in the FAM and another individual follows."<<endl;
  cerr<<"\t                 This list is also the pool for --nsample. Without --nsample, LD scores"<<endl;
  cerr<<"\t                 use the same individuals. Default: every individual in the FAM."<<endl;
  cerr<<"\t--nsample      : Number of individuals drawn without replacement for LD scores."<<endl;
  cerr<<"\t                 The draw is from --keep, or from the full FAM when --keep is absent."<<endl;
  cerr<<"\t                 The same individuals are used on every chromosome. If the pool has this"<<endl;
  cerr<<"\t                 many people or fewer, everyone in the pool is used. Allele frequency and"<<endl;
  cerr<<"\t                 missingness still use the full pool."<<endl;
  cerr<<"\t--write-sample : Write the LD-score individuals to [prefix].sample. Two columns, FID and IID,"<<endl;
  cerr<<"\t                 with a header. The same list --keep reads. With --nsample this is the draw,"<<endl;
  cerr<<"\t                 and without it this is the eligible set."<<endl;
  cerr<<"\t--seed         : Seed for --nsample. An integer from 0 through 4294967295. Default 1."<<endl;
  cerr<<"\t                 Used only when --nsample draws a subset."<<endl;
  cerr<<"\t--nthread      : POSIX threads for the LD-score SNP loop within each chromosome. Default 1."<<endl;
  cerr<<"\t                 Threads share the chromosome genotype matrix. One thread keeps the"<<endl;
  cerr<<"\t                 single-thread scores. More threads can change the last bits."<<endl;
  cerr<<"\t--map-average  : One sex-averaged genetic map for autosomes 1-22. Five columns, in order:"<<endl;
  cerr<<"\t                 chromosome, begin bp, end bp, cM per Mb, cumulative cM at end bp."<<endl;
  cerr<<"\t                 A header is optional and is not matched by name."<<endl;
  cerr<<"\t                 Chromosomes in both the map and the BIM are scored, in order from 1 to 22."<<endl;
  cerr<<"\t                 X, Y, MT and any other chromosome are dropped. BIM column 3 (cM) is ignored."<<endl;
  cerr<<"\t                 Not used with --mb-window."<<endl;
  cerr<<"\t--cm-window    : LD window in centiMorgans. Pairs at or beyond this distance are excluded."<<endl;
  cerr<<"\t                 Default 1. Not used with --mb-window."<<endl;
  cerr<<"\t--mb-window    : LD window in megabases. Pairs at or beyond this distance are excluded."<<endl;
  cerr<<"\t                 --mb-window 1 keeps pairs less than 1,000,000 bp apart."<<endl;
  cerr<<"\t                 --map-average, --cm-window, and --interpolate are not used."<<endl;
  cerr<<"\t                 Every autosome in the BIM is scored."<<endl;
  cerr<<"\t--extract      : Text file of SNP identifiers, one per line. Only BIM SNPs with a"<<endl;
  cerr<<"\t                 listed identifier are analysed. Blank lines and lines starting with #"<<endl;
  cerr<<"\t                 are ignored. A repeated identifier is kept once. The cut is applied"<<endl;
  cerr<<"\t                 before missingness, MAF, and the map or base-pair window."<<endl;
  cerr<<"\t                 Use either --extract or --mextract."<<endl;
  cerr<<"\t--mextract     : Text file of SNP-list paths, one per line. Each path is an --extract"<<endl;
  cerr<<"\t                 file. Identifiers from every file are pooled. Blank lines and lines"<<endl;
  cerr<<"\t                 starting with # are ignored. A repeated path is an error."<<endl;
  cerr<<"\t--maf          : Minimum minor allele frequency (fraction). SNPs below this value are dropped."<<endl;
  cerr<<"\t                 Default 0.01."<<endl;
  cerr<<"\t--max-missing  : Maximum missingness (fraction, e.g. 0.05 = 5%). Default 0.05. Must be in [0,1]."<<endl;
  cerr<<"\t--interpolate  : Place each SNP by linear interpolation inside its map interval. Off by default."<<endl;
  cerr<<"\t                 Off: a SNP inside an interval takes the cM at Begin, and a SNP outside the span is dropped."<<endl;
  cerr<<"\t                 On: interior SNPs are interpolated, and SNPs outside the span are clamped to the end values."<<endl;
  cerr<<"\t--out          : Prefix for [prefix].ldscore.txt and [prefix].log. Default is ldscore."<<endl;
  cerr<<"\t                 LD scores have four columns: chromosome, SNP id, base pair, LD score."<<endl;
}

} // namespace

int main(int argc, char *argv[]){
  string bfile;
  string mbfile;
  string keepFile;
  string extract_file;
  string mextract_file;
  string mapFile;
  string outPrefix = "ldscore";
  bool do_interpolate = false;
  double window_cm = 1.0;
  double mb_window = 0.0;
  bool mb_set = false;
  double min_maf = 0.01;
  double max_missing = 0.05;
  uint32_t nsample = 0;
  uint32_t seed = 1;
  int nthread = 1;
  bool nsample_set = false;
  bool seed_set = false;
  bool write_sample = false;

  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    return 1;
  }
  if(string(argv[1])=="--help"){
    print_help();
    return 1;
  }

  for(int i=1;i<argc;++i){
    const string sw = argv[i];
    auto need = [&](const string &name){
      if(i+1>=argc){
        cerr<<"\tMissing value for "<<name<<"."<<endl;
        exit(1);
      }
      return string(argv[++i]);
    };
    if(sw=="--bfile") bfile = need(sw);
    else if(sw=="--mbfile") mbfile = need(sw);
    else if(sw=="--keep") keepFile = need(sw);
    else if(sw=="--extract") extract_file = need(sw);
    else if(sw=="--mextract") mextract_file = need(sw);
    else if(sw=="--nsample"){
      const string tok = need(sw);
      if(!parse_full_u32(tok, nsample) || nsample<1){
        cerr<<"\t--nsample must be a positive integer."<<endl;
        return 1;
      }
      nsample_set = true;
    }
    else if(sw=="--seed"){
      const string tok = need(sw);
      if(!parse_full_u32(tok, seed)){
        cerr<<"\t--seed must be an integer from 0 through 4294967295."<<endl;
        return 1;
      }
      seed_set = true;
    }
    else if(sw=="--nthread"){
      const string tok = need(sw);
      uint32_t nt = 0;
      if(!parse_full_u32(tok, nt) || nt<1 || nt>static_cast<uint32_t>(numeric_limits<int>::max())){
        cerr<<"\t--nthread must be an integer >= 1."<<endl;
        return 1;
      }
      nthread = static_cast<int>(nt);
    }
    else if(sw=="--write-sample") write_sample = true;
    else if(sw=="--map-average") mapFile = need(sw);
    else if(sw=="--map"){
      cerr<<"\t--map has been renamed. Use --map-average."<<endl;
      return 1;
    }
    else if(sw=="--cm-window"){
      const string tok = need(sw);
      if(!parse_full_double(tok, window_cm)){
        cerr<<"\t--cm-window must be a positive number of centiMorgans."<<endl;
        return 1;
      }
    }
    else if(sw=="--mb-window"){
      const string tok = need(sw);
      if(!parse_full_double(tok, mb_window) || !(mb_window>0.0)){
        cerr<<"\t--mb-window must be a positive number of megabases."<<endl;
        return 1;
      }
      mb_set = true;
    }
    else if(sw=="--maf"){
      const string tok = need(sw);
      if(!parse_full_double(tok, min_maf)){
        cerr<<"\t--maf must be a fraction between 0 and 0.5."<<endl;
        return 1;
      }
    }
    else if(sw=="--max-missing"){
      const string tok = need(sw);
      if(!parse_full_double(tok, max_missing)){
        cerr<<"\t--max-missing must be a fraction in [0, 1] (e.g. 0.05 for 5%)."<<endl;
        return 1;
      }
    }
    else if(sw=="--interpolate") do_interpolate = true;
    else if(sw=="--out") outPrefix = need(sw);
    else if(sw=="--help"){
      print_help();
      return 1;
    }
    else{
      cerr<<"\tUnknown argument "<<sw<<". Type --help for more details."<<endl;
      return 1;
    }
  }

  if(!bfile.empty() && !mbfile.empty()){
    cerr<<"\tUse either --bfile or --mbfile, not both."<<endl;
    return 1;
  }
  if(bfile.empty() && mbfile.empty()){
    cerr<<"\tA PLINK prefix must be specified. [Use: --bfile or --mbfile]"<<endl;
    return 1;
  }
  if(!extract_file.empty() && !mextract_file.empty()){
    cerr<<"\tUse either --extract or --mextract, not both."<<endl;
    return 1;
  }
  if(!mb_set && mapFile.empty()){
    cerr<<"\tA genetic map must be specified. [Use: --map-average]"<<endl;
    return 1;
  }
  if(!mb_set && !(window_cm>0.0)){
    cerr<<"\t--cm-window must be a positive number of centiMorgans."<<endl;
    return 1;
  }
  if(!(min_maf>=0.0 && min_maf<=0.5)){
    cerr<<"\t--maf must be a fraction between 0 and 0.5."<<endl;
    return 1;
  }
  if(!(max_missing>=0.0 && max_missing<=1.0)){
    cerr<<"\t--max-missing must be a fraction in [0, 1] (e.g. 0.05 for 5%)."<<endl;
    return 1;
  }

  vector<string> prefixes;
  if(!mbfile.empty()){
    string err;
    if(!read_mbfile(mbfile, prefixes, err)){
      cerr<<err<<endl;
      return 1;
    }
  }else{
    prefixes.push_back(bfile);
  }
  const char *exts[] = {".bed", ".bim", ".fam"};
  for(const string &pref : prefixes){
    for(const char *ext : exts){
      const string path = pref + ext;
      ifstream in(path.c_str());
      if(!in){
        cerr<<"Error reading "<<path<<endl;
        return 1;
      }
    }
  }

  unordered_set<string> extract_ids;
  int n_extract_lines = 0;
  int n_extract_files = 0;
  const bool use_extract = !extract_file.empty() || !mextract_file.empty();
  if(use_extract){
    string err;
    vector<string> lists;
    if(!extract_file.empty()){
      lists.push_back(extract_file);
    }else if(!read_extract_list(mextract_file, lists, err)){
      cerr<<err<<endl;
      return 1;
    }
    n_extract_files = static_cast<int>(lists.size());
    for(const string &path : lists){
      if(!read_snp_list(path, extract_ids, n_extract_lines, err)){
        cerr<<err<<endl;
        return 1;
      }
    }
    if(extract_ids.empty()){
      cerr<<"No SNP identifiers in "<<(extract_file.empty() ? mextract_file : extract_file)<<"."<<endl;
      return 1;
    }
  }

  const string famfile = prefixes.front() + ".fam";
  const string logFile = outPrefix + ".log";
  const string ldFile  = outPrefix + ".ldscore.txt";
  const string sampleFile = outPrefix + ".sample";

  const auto t0 = chrono::steady_clock::now();
  time_t t = time(0);
  struct tm *now = localtime(&t);

  ofstream fileLog(logFile.c_str());
  if(!fileLog){
    cerr<<"Error writing log file "<<logFile<<endl;
    return 1;
  }
  SayFn say;
  say.log = &fileLog;

  vector<string> fam_fid, fam_iid;
  int n_dup_iid = 0;
  {
    string err;
    if(!load_fam(famfile, fam_fid, fam_iid, n_dup_iid, err)){
      cerr<<err<<endl;
      return 1;
    }
  }
  const int n = static_cast<int>(fam_iid.size());
  for(size_t pi=1;pi<prefixes.size();++pi){
    vector<string> fid2, iid2;
    int dup2 = 0;
    string err;
    const string fam2 = prefixes[pi] + ".fam";
    if(!load_fam(fam2, fid2, iid2, dup2, err)){
      cerr<<err<<endl;
      return 1;
    }
    if(fid2!=fam_fid || iid2!=fam_iid){
      cerr<<"FAM "<<fam2<<" does not match "<<famfile<<". Individuals must be the same and in the same order."<<endl;
      return 1;
    }
  }
  vector<int> sample_ids;
  vector<string> keep_notes;
  {
    string err;
    if(!select_samples(fam_fid, fam_iid, keepFile, sample_ids, keep_notes, err)){
      cerr<<err<<endl;
      return 1;
    }
  }
  bool ld_drawn = false;
  vector<int> ld_ids = sample_ids;
  if(nsample_set) ld_ids = draw_ld_sample(sample_ids, nsample, seed, ld_drawn);
  if(write_sample){
    ofstream samp(sampleFile.c_str());
    if(!samp){
      cerr<<"Error writing "<<sampleFile<<endl;
      return 1;
    }
    samp<<"FID\tIID\n";
    for(int idx : ld_ids) samp<<fam_fid[idx]<<"\t"<<fam_iid[idx]<<"\n";
    if(!samp){
      cerr<<"Error writing "<<sampleFile<<endl;
      return 1;
    }
  }

  const int numBytes = (n + PACK_DENSITY - 1) / PACK_DENSITY;
  vector<BimSet> bims;
  bims.reserve(prefixes.size());
  for(const string &pref : prefixes){
    BimSet one;
    one.prefix = pref;
    one.bedfile = pref + ".bed";
    const string bimfile = pref + ".bim";
    string err;
    if(!load_bim(bimfile, one.snps, err)){
      cerr<<err<<endl;
      return 1;
    }
    if(!index_bim_autosomes(bimfile, one.snps, one, err)){
      cerr<<err<<endl;
      return 1;
    }
    if(!check_bed_size(one.bedfile, one.snps.size(), numBytes, err)){
      cerr<<err<<endl;
      return 1;
    }
    bims.push_back(std::move(one));
  }
  unordered_set<string> extract_empty;
  ExtractCounts extract_counts;
  if(use_extract){
    apply_extract(bims, extract_ids, extract_counts, extract_empty);
  }

  unordered_map<string, vector<MapPoint>> map_auto;
  vector<string> skipped_map;
  if(!mb_set){
    string err;
    if(!read_interval_map(mapFile, map_auto, skipped_map, err)){
      cerr<<err<<endl;
      return 1;
    }
  }
  const double window_pos = mb_set ? mb_window * 1e6 : window_cm;

  ofstream ldout(ldFile.c_str());
  if(!ldout){
    cerr<<"Error writing "<<ldFile<<endl;
    return 1;
  }
  ldout<<"CHR\tSNP\tBP\tL2\n";

  {
    ostringstream oss;
    oss<<">>> LD score calculation <<<\n";
    oss<<"# Analysis starts : "<<(now->tm_year+1900)<<'-'<<(now->tm_mon+1)<<'-'<<now->tm_mday
       <<" at "<<now->tm_hour<<":"<<now->tm_min<<":"<<now->tm_sec<<".\n";
    oss<<"# Command:";
    for(int i=0;i<argc;++i) oss<<" "<<argv[i];
    oss<<"\n";
    if(prefixes.size()==1){
      oss<<"# Genotype file: [BED: "<<prefixes[0]<<".bed], [BIM: "<<prefixes[0]<<".bim], [FAM: "<<prefixes[0]<<".fam]\n";
    }else{
      oss<<"# Genotype list: "<<mbfile<<" ("<<prefixes.size()<<" prefixes)\n";
      for(const string &pref : prefixes){
        oss<<"#   "<<pref<<".bed/.bim/.fam\n";
      }
    }
    if(!keepFile.empty()) oss<<"# Keep file: "<<keepFile<<"\n";
    if(use_extract){
      oss<<"# SNP list"<<((n_extract_files==1) ? "" : "s")<<": "
         <<(extract_file.empty() ? mextract_file : extract_file)
         <<" ("<<n_extract_files<<" file"<<(n_extract_files==1 ? "" : "s")<<", "
         <<extract_ids.size()<<" identifier"<<(extract_ids.size()==1 ? "" : "s");
      if(n_extract_lines != static_cast<int>(extract_ids.size())){
        oss<<" from "<<n_extract_lines<<" lines";
      }
      oss<<").\n";
    }
    if(write_sample) oss<<"# LD sample: "<<sampleFile<<"\n";
    if(mb_set){
      ostringstream bpss;
      if(window_pos == floor(window_pos) && window_pos < 1e15) bpss<<static_cast<unsigned long long>(llround(window_pos));
      else bpss<<window_pos;
      oss<<"# LD window: "<<mb_window<<" Mb ("<<bpss.str()<<" bp).\n";
      oss<<"# --map-average and --cm-window are not used.\n";
      if(do_interpolate) oss<<"# --interpolate is not used with --mb-window.\n";
    }else{
      oss<<"# Average map: "<<mapFile<<"\n";
    }
    oss<<"# Output LD scores: "<<ldFile<<"\n";
    oss<<"# cm-window="<<window_cm
       <<" maf="<<min_maf
       <<" max-missing="<<max_missing
       <<" interpolate="<<(do_interpolate ? "on" : "off")
       <<" nthread="<<nthread;
    if(mb_set) oss<<" mb-window="<<mb_window;
    if(nsample_set) oss<<" nsample="<<nsample<<" seed="<<seed;
    say(oss.str());
  }
  if(use_extract){
    say("# Extract: kept "+to_string(extract_counts.n_kept)+" autosomal BIM SNP"+(extract_counts.n_kept==1 ? "" : "s")+
        ", dropped "+to_string(extract_counts.n_dropped)+"; "+
        to_string(extract_counts.n_absent)+" listed identifier"+(extract_counts.n_absent==1 ? " is" : "s are")+
        " absent from the BIM.");
    for(int c=1;c<=22;++c){
      const string key = to_string(c);
      if(extract_empty.count(key)){
        say("# [Note] Chromosome "+key+" has no SNPs in the SNP list; skipping.");
      }
    }
    if(extract_counts.n_kept==0){
      cerr<<"No autosomal BIM SNPs remain after "<<(extract_file.empty() ? "--mextract" : "--extract")<<"."<<endl;
      return 1;
    }
  }
  say("# Found "+to_string(n)+" individuals in "+famfile+
      "; "+to_string(sample_ids.size())+" eligible for allele frequency and missingness.");
  if(prefixes.size()>1){
    say("# Checked "+to_string(prefixes.size())+" FAM files: same individuals in the same order.");
  }
  for(const string &note : keep_notes) say(note);
  if(nsample_set && ld_drawn){
    say("# --nsample "+to_string(nsample)+" seed "+to_string(seed)+
        ": drew "+to_string(ld_ids.size())+" of "+to_string(sample_ids.size())+
        " eligible individuals for LD scores.");
  }else if(nsample_set){
    say("# --nsample "+to_string(nsample)+": "+to_string(sample_ids.size())+
        " eligible individuals, using all of them for LD scores. Seed was not applied.");
  }else{
    say("# LD scores use all "+to_string(sample_ids.size())+" eligible individuals.");
    if(seed_set) say("# [Note] --seed is unused without --nsample.");
  }
  if(write_sample){
    say("# Wrote "+to_string(ld_ids.size())+" LD-score individuals to "+sampleFile+".");
  }
  if(ld_ids.size()<3){
    say("# [Note] Fewer than 3 individuals are used for LD scores; pairwise r^2 is 0 and LD scores stay at 1.");
  }
  if(n_dup_iid>0){
    say("# [Warning] "+to_string(n_dup_iid)+" duplicate IID(s) in "+famfile+".");
  }
  if(!mb_set && !skipped_map.empty()){
    say("# [Note] Dropped non-autosomal map chromosomes: "+join_keys(skipped_map)+".");
  }
  {
    vector<string> dropped;
    unordered_set<string> seen_drop;
    for(const BimSet &one : bims){
      for(const string &key : one.dropped){
        if(seen_drop.insert(key).second) dropped.push_back(key);
      }
    }
    if(!dropped.empty()){
      say("# [Note] Dropped non-autosomal BIM chromosomes: "+join_keys(dropped)+".");
    }
  }

  unordered_set<string> bim_chrs;
  for(const BimSet &one : bims){
    for(const auto &kv : one.chr_snps) bim_chrs.insert(kv.first);
  }
  vector<string> map_only, bim_only;
  if(!mb_set){
    for(int c=1;c<=22;++c){
      const string key = to_string(c);
      const bool in_map = map_auto.count(key)>0;
      const bool in_bim = bim_chrs.count(key)>0;
      if(in_map && !in_bim) map_only.push_back(key);
      if(in_bim && !in_map) bim_only.push_back(key);
    }
  }
  if(!map_only.empty()){
    say("# [Warning] Map chromosomes with no BIM SNPs: "+join_keys(map_only)+".");
  }
  if(!bim_only.empty()){
    say("# [Warning] BIM chromosomes missing from the map: "+join_keys(bim_only)+".");
  }

  struct ChrJob { int bim_i; string key; };
  vector<ChrJob> jobs;
  unordered_map<string, vector<string>> chunk_from;
  for(int bi=0; bi<static_cast<int>(bims.size()); ++bi){
    for(int c=1;c<=22;++c){
      const string key = to_string(c);
      if(!bims[bi].chr_snps.count(key)) continue;
      if(!mb_set && !map_auto.count(key)) continue;
      const string where = prefixes.size()>1 ? " in "+bims[bi].prefix : "";
      if(mb_set){
        say("# Chromosome "+bims[bi].chr_label[key]+where+": "+to_string(bims[bi].chr_snps[key].size())+" BIM SNPs.");
      }else{
        say("# Chromosome "+bims[bi].chr_label[key]+where+": "+to_string(bims[bi].chr_snps[key].size())+
            " BIM SNPs, "+to_string(map_auto[key].size())+" map boundaries.");
      }
      chunk_from[key].push_back(bims[bi].prefix);
      jobs.push_back(ChrJob{bi, key});
    }
  }
  for(int c=1;c<=22;++c){
    const string key = to_string(c);
    if(chunk_from[key].size()<2) continue;
    say("# [Note] Chromosome "+key+" is in "+to_string(chunk_from[key].size())+
        " prefixes ("+join_keys(chunk_from[key])+"). Each chunk is scored separately.");
  }
  if(jobs.empty()){
    cerr<<(mb_set ? "No autosomes in the BIM." : "No autosomes overlap the map and the BIM.")<<endl;
    return 1;
  }

  vector<char> packed(numBytes);
  vector<char> unpacked(static_cast<size_t>(numBytes) * PACK_DENSITY);
  int n_written = 0;
  int n_floor = 0;
  unordered_set<string> wrote_chr;
  unordered_map<string, unordered_set<string>> written_ids;

  for(const ChrJob &job : jobs){
    const BimSet &bim = bims[job.bim_i];
    const string &key = job.key;
    const string label = prefixes.size()>1 ? bim.chr_label.at(key)+" in "+bim.prefix : bim.chr_label.at(key);
    say("# Dosage matrix upper bound for chromosome "+label+": "+
        mem_line(ld_ids.size(), bim.chr_snps.at(key).size(), 1)+
        (mb_set ? " before MAF and missingness filters."
                : " before MAF, missingness, and map filters."));
    ChrResult scored;
    string err;
    const vector<MapPoint> *gmap = mb_set ? nullptr : &map_auto.at(key);
    if(!score_chromosome(bim.bedfile, bim.snps, bim.chr_snps.at(key), gmap,
                         sample_ids, ld_ids,
                         do_interpolate, min_maf, max_missing,
                         numBytes, packed.data(), unpacked.data(), scored, err)){
      cerr<<err<<endl;
      return 1;
    }
    say("# SNP filters: dropped "+to_string(scored.n_maf)+" for MAF, "+
        to_string(scored.n_missed)+" for missingness"+
        (mb_set || do_interpolate ? string("") : (", "+to_string(scored.n_nomap)+" outside the map"))+
        "; kept "+to_string(scored.snps.size())+".");
    if(scored.snps.empty()){
      say("# [Warning] No SNPs remaining after filters on chromosome "+label+"; skipping.");
      continue;
    }
    if(mb_set){
      ostringstream oss;
      oss<<"# Kept SNP bp range: "<<scored.snps.front().bp<<" to "<<scored.snps.back().bp
         <<" (span "<<(scored.snps.back().bp-scored.snps.front().bp)<<" bp).";
      say(oss.str());
    }else{
      const char *mid = do_interpolate ? " interpolated=" : " left-edge=";
      ostringstream oss;
      oss<<"# Map onto BIM bp: boundary="<<scored.n_exact<<mid<<scored.n_interp
         <<" clamped="<<scored.n_clamp<<".\n"
         <<"# Kept SNP cM range: "<<scored.cm.front()<<" to "<<scored.cm.back()
         <<" (span "<<(scored.cm.back()-scored.cm.front())<<" cM); map bp "
         <<map_auto.at(key).front().bp<<" to "<<map_auto.at(key).back().bp<<".";
      say(oss.str());
    }
    say("# LD genotypes: "+mem_line(ld_ids.size(), scored.snps.size(), 1)+".");
    const int n_sample = static_cast<int>(ld_ids.size());
    const int M = static_cast<int>(scored.snps.size());
    int n_used = 1;
    if(nthread>1 && M>1){
      const int nblocks = (M + LD_SNP_BLOCK - 1) / LD_SNP_BLOCK;
      n_used = nthread < nblocks ? nthread : nblocks;
    }
    {
      ostringstream win;
      if(mb_set) win<<mb_window<<" Mb";
      else win<<window_cm<<" cM";
      say("# Computing LD scores on "+to_string(n_sample)+
          " individuals ("+win.str()+" window) using "+to_string(n_used)+" thread"+
          (n_used==1 ? "" : "s")+".");
    }
    vector<double> ldscores;
    {
      string err;
      if(!compute_ld_scores_threaded(scored.gts.data(), n_sample, scored.cm, window_pos,
                                     nthread, ldscores, n_used, err)){
        cerr<<err<<endl;
        return 1;
      }
    }
    scored.gts.clear();
    scored.gts.shrink_to_fit();
    scored.l2.resize(M);
    for(int k=0;k<M;++k){
      if(!(ldscores[k] >= L2_FLOOR)){
        ldscores[k] = L2_FLOOR;
        scored.n_floor++;
      }
      scored.l2[k] = ldscores[k];
    }
    double sum = 0.0;
    for(int k=0;k<static_cast<int>(scored.snps.size());++k){
      if(!written_ids[key].insert(scored.snps[k].id).second){
        cerr<<"Duplicate SNP id \""<<scored.snps[k].id<<"\" on chromosome "<<key<<" across prefixes."<<endl;
        return 1;
      }
      ldout<<scored.snps[k].chr<<"\t"<<scored.snps[k].id<<"\t"
           <<scored.snps[k].bp<<"\t"<<scored.l2[k]<<"\n";
      sum += scored.l2[k];
    }
    if(!ldout){
      cerr<<"Error writing "<<ldFile<<endl;
      return 1;
    }
    n_written += static_cast<int>(scored.snps.size());
    n_floor += scored.n_floor;
    wrote_chr.insert(key);
    {
      ostringstream mean;
      mean<<(sum / static_cast<double>(scored.snps.size()));
      say("# Chromosome "+label+": wrote "+to_string(scored.snps.size())+
          " LD scores (mean "+mean.str()+").");
    }
  }

  if(n_written==0){
    cerr<<"No SNPs remaining after filters."<<endl;
    return 1;
  }
  if(n_floor>0){
    say("# Floored "+to_string(n_floor)+" LD score"+(n_floor==1 ? "" : "s")+
        " at "+to_string(L2_FLOOR)+".");
  }

  const double elapsed = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
  time_t t2 = time(0);
  struct tm *now2 = localtime(&t2);
  {
    ostringstream oss;
    const int n_chr_written = static_cast<int>(wrote_chr.size());
    oss<<"# Wrote "<<n_written<<" LD scores across "<<n_chr_written
       <<" autosome"<<(n_chr_written==1 ? "" : "s")<<" to "<<ldFile<<".\n";
    oss<<"# Analysis ends: "<<(now2->tm_year+1900)<<'-'<<(now2->tm_mon+1)<<'-'<<now2->tm_mday
       <<" at "<<now2->tm_hour<<":"<<now2->tm_min<<":"<<now2->tm_sec<<".\n";
    oss<<"# Time elapsed: "<<elapsed<<" seconds.\n";
    oss<<"<<< LD score calculation >>>";
    say(oss.str());
  }
  return 0;
}
