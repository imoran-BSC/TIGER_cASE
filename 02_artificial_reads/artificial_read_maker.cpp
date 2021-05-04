/*
 *  ARTIFICIAL_READ_MAKER.CPP
 *
 *  Reads a list of SNPs and a fasta genome, and generates
 *  all possible reads overlapping the SNPs, a histogram of
 *  the SNP density in each base, and a list of the discarded
 *  excessively polymorphic regions
 */

#include "../99_utils/class_fstream.h"


using namespace std;
const bool debug = 0;  // DEBUG FLAG
const int line_char_count = 50;  // Char width of .fasta files
const double min_maf_value = 0.005;  // MAF threshold for exclusion
double acceptable_snp_ratio;  // Max SNP ratio to generate all the possible reads
unsigned read_length, max_num_snps, read_counter=0;

string funct_to_uppercase( string );  // Needs to be here because some classes use it

/*
 *   CLASSES
 */

class class_gtf
{  // Format GTF - used to get all gene and lncRNA exons and their junctions
public:
  unsigned ini, fin, length;
  string chr, type, name, key;

  set<string> snp_set;

  class_gtf() : ini(0), fin(0), length(0),
    chr(""), type(""), name(""), key("") {}

  void get_gtf( class_fstream& in )
    {
    string dummy;
    stringstream ss;

    in.file >> chr >> dummy >> type >> ini >> fin
      >> dummy >> dummy >> dummy >> name;

    length = fin-ini+1;
    ss << chr << ":" << ini << "-" << fin;
    key = ss.str();
    }
};

class class_snp
{  // Info of a dbSNP in 1kG format
public:
  unsigned pose;
  double maf;
  char ref, alt;
  string chr, name, sref, salt;

  class_snp() : pose(0), ref('\0'), alt('\0'),
    chr(""), name(""), sref(""), salt("") {}

  void get_snp( class_fstream& in )
    {
    size_t p0, p1, p2;
    string nchr, dummy;

    in.file >> nchr >> name
      >> dummy >> dummy >> maf >> dummy;

    p0 = name.find(":");
    p1 = name.find("_");
    p2 = name.find("_", p1+1);

    chr = "chr" + nchr;
    sref = name.substr( p1+1, p2-p1-1 );
    salt = name.substr( p2+1 );
    pose = atol( name.substr( p0+1, p1-p0-1 ).c_str() );
    }

  bool isgood()
    {  // Check if the SNP has a proper maf and its binary
    return maf >= min_maf_value && sref.length() == 1 && salt.length() == 1;
    }

  void make_nt_chars()
    {  // Transform the string ref/alt nts to chars - REQUIRES isgood()!
    ref = sref[0]; alt = salt[0];
    sref = ""; salt = "";  // and remove the strings from mem
    }
};

class class_newread
{  // Methods to create, update and write new readSNPs
public:
  unsigned status, ini_pose;
  string chr, nucleotides;

  vector<class_snp> snp_vec;
  map<string, class_gtf> exons_map;

  class_newread() : status(0), ini_pose(0), chr(""), nucleotides("") {}

  void create( class_snp& s, unsigned& lcount, string& str )
    {  // Creates a newread with the start info and firts overlapping SNP
    ini_pose = s.pose - read_length; chr = s.chr;
    nucleotides = str.substr( ini_pose - (lcount-1)*line_char_count );
    snp_vec.clear(); snp_vec.push_back(s);
    }

  void update( string& str )
    {  // Update the nucleotide string
    int last_snp_pose = snp_vec.back().pose - ini_pose;

    nucleotides = nucleotides + str;
    if( (int)nucleotides.size() - last_snp_pose > (int)read_length )
      { nucleotides = funct_to_uppercase( nucleotides.substr(0,last_snp_pose + read_length) );
        status = 1; } // Poised for writting
    }

  void write_reads( class_fstream& out, class_fstream& forbidden, map<int, int>& snp_density_map )
    {  // With the nucleotide chain and all overlapping SNPs,
       // creates all possible reads and writes them out, plus
       // makes an hustogram of the number of SNPs / read_length.
    bool forbidden_on=0, contained_in_exon=0;
    unsigned vec_indx=0, nt_vec_current_size, asnps_size,
      vec_len = snp_vec.size(), str_length = nucleotides.size();
    string nt_tmp, read_tmp;
    stringstream ss;

    class_gtf forbidden_reg;

    vector<string> nt_vec;
    list<class_snp> active_snps;
    list<class_snp>::iterator asnps_it;
    set<string>::iterator  set_it;
    map<int,int>::iterator map_it;
    map<string, class_gtf>::iterator exon_it;

    for(unsigned k=0; k<(str_length-read_length); k++)
      {
      contained_in_exon=0;
      for( asnps_it = active_snps.begin(); asnps_it != active_snps.end(); )  // Remove past SNPs
        { if( (*asnps_it).pose <= ini_pose+k ) asnps_it = active_snps.erase( asnps_it );
          else ++asnps_it; }

      for(exon_it = exons_map.begin(); exon_it != exons_map.end(); ++exon_it )
        if( ini_pose+k >= (*exon_it).second.ini &&
            ini_pose+k+read_length <= (*exon_it).second.fin )
          { contained_in_exon = 1; break; } // Read completely contained within exon

      while( vec_indx < vec_len &&
             snp_vec[vec_indx].pose <= k+ini_pose+read_length )  // Add active SNPs
        { if(forbidden_on) forbidden_reg.snp_set.insert(snp_vec[vec_indx].name);
          active_snps.push_back(snp_vec[vec_indx]); vec_indx++; }
      asnps_size = active_snps.size();

      map_it = snp_density_map.find( asnps_size );  // Histogram of SNP density
      if( map_it == snp_density_map.end() ) snp_density_map[ asnps_size ] = 1;
      else (*map_it).second++;

      if( contained_in_exon )
        {
        if( asnps_size <= read_length*acceptable_snp_ratio )
          {  // Process reads if less than XX% of the read is made of SNPs
          contained_in_exon = 0;
          if(forbidden_on)
            {  // End of an open forbidden region
            long int far_exon_border=0;
            forbidden_on = 0;

            for(exon_it = exons_map.begin(); exon_it != exons_map.end(); ++exon_it )
              if( (*exon_it).second.fin > far_exon_border )
                far_exon_border = (*exon_it).second.fin;
            if( ini_pose+k+read_length-1 > far_exon_border ) forbidden_reg.fin = far_exon_border;
            else forbidden_reg.fin = ini_pose+k+read_length-1;

            forbidden.file << forbidden_reg.chr << "\t" << forbidden_reg.ini << "\t"
              << forbidden_reg.fin << "\t" << forbidden_reg.snp_set.size() << "\t";
            for( set_it=forbidden_reg.snp_set.begin(); set_it!=forbidden_reg.snp_set.end(); ++set_it )
              forbidden.file << (*set_it) << ",";
            forbidden.file << "\n";
            forbidden_reg.snp_set.clear();
            }

          nt_vec.clear();
          nt_vec.push_back( nucleotides.substr(k,read_length) );  // This is the reference read

          for( asnps_it = active_snps.begin(); asnps_it != active_snps.end(); ++asnps_it )
            {  // Create all alternative reads, looping through the active snps in the fragment
            nt_vec_current_size = nt_vec.size();
            for(unsigned j=0; j<nt_vec_current_size; j++)
              {
              nt_tmp = nt_vec[j];
              // Sanity check
              if( nt_tmp[ (*asnps_it).pose-ini_pose-1-k ] != (*asnps_it).ref )
                { cout << "We got a mismatch w the reference genome:\n "
                    << nt_tmp << "\n with " << (*asnps_it).ref
                    << " at pose " << (*asnps_it).pose-ini_pose-1-k
                    << "(" << (*asnps_it).pose << ")\n"; cin.get(); }
              nt_tmp[ (*asnps_it).pose-ini_pose-1-k ] = (*asnps_it).alt;
              nt_vec.push_back( nt_tmp );
              }
            }

          nt_vec_current_size = nt_vec.size();
          for(unsigned j=0; j<nt_vec_current_size; j++)
            {  // Write reads to output in fasta phred33 format
            ss.str(""); read_tmp = "";
            // Read header, contains info about where was it generated,
            // number of overlapping SNPs and of the read itself
            ss << "@" << chr << ":" << ini_pose+k+1 << ":" << pow(2,asnps_size) << ":" << j;
            read_tmp = ss.str() + "\n" + nt_vec[j] + "\n+\n";
            for(unsigned i=0; i<read_length; i++) read_tmp += "I";  // Quality, in phred33 its almost max

            out.file << read_tmp << "\n"; read_counter++;
            }
          }
        else if( !forbidden_on )
          {
          forbidden_on=1;
          forbidden_reg.chr = chr; forbidden_reg.ini = ini_pose+k+1;
          for( asnps_it = active_snps.begin(); asnps_it != active_snps.end(); ++asnps_it )
            forbidden_reg.snp_set.insert( (*asnps_it).name );  // List of SNPs in the region
          }
        }
      }

    if(forbidden_on)
      {  // End of an open forbidden region due to end of region
      long int far_exon_border=0;
      forbidden_on = 0;

      for(exon_it = exons_map.begin(); exon_it != exons_map.end(); ++exon_it )
        if( (*exon_it).second.fin > far_exon_border )
          far_exon_border = (*exon_it).second.fin;
      if( ini_pose+str_length > far_exon_border ) forbidden_reg.fin = far_exon_border;
      else forbidden_reg.fin = ini_pose+str_length;

      forbidden.file << forbidden_reg.chr << "\t" << forbidden_reg.ini << "\t"
        << forbidden_reg.fin << "\t";
      for( set_it=forbidden_reg.snp_set.begin(); set_it!=forbidden_reg.snp_set.end(); ++set_it )
        forbidden.file << (*set_it) << ",";
      forbidden.file << "\t" << forbidden_reg.snp_set.size() << "\n";
      forbidden_reg.snp_set.clear();
      }
    }
};

class class_junction
{  // Junction info
public:
  bool active_leftmost, close_left, close_right;
  unsigned leftmost, rightmost;  // Last and first exon bases
  string chr, left_nt, right_nt;

  class_gtf left_exon, right_exon;

  vector<class_snp> snp_vec;
  map<string, class_gtf> left_exons_map, right_exons_map;

  class_junction() : active_leftmost(0), close_left(0), close_right(0), leftmost(0), rightmost(0),
    chr(""), left_nt(""), right_nt("") {}

  void make_junction( class_gtf& gtf_left, class_gtf& gtf_right )
    { chr = gtf_left.chr; leftmost = gtf_left.fin; rightmost = gtf_right.ini; }

  void update_left_nt( string& str, unsigned& lcount )
    {  // Init or update left nt_string
    if(left_nt == "")
      left_nt = str.substr( leftmost-read_length+1 - (lcount-1)*line_char_count, string::npos );
    else left_nt = left_nt + str;
    }

  void close_left_nt( string& str )
    {  // Finalizes left nt
    left_nt = left_nt + str;
    left_nt = funct_to_uppercase( left_nt.substr(0, read_length-1) );
    }

  void update_right_nt( string& str, unsigned& lcount )
    {  // Init or update right nt_string
    if(right_nt == "")
      right_nt = str.substr( rightmost - (lcount-1)*line_char_count-1, string::npos );
    else right_nt = right_nt + str;
    }

  void close_right_nt( string& str )
    {  // Finalizes right nt
    right_nt = right_nt + str;
    right_nt = funct_to_uppercase( right_nt.substr(0, read_length-1) );
    }

  void write( class_fstream& out )
    {  // Create all permutations and writeout
    int exon_limit=0;  // 0=good, 1=left, 2=right limit, 3=both, 4=both < rlen
    unsigned vec_indx=0, nt_vec_current_size, asnps_size,
      vec_len = snp_vec.size(), str_length, ini_pose = leftmost-read_length+1;
    string nt_tmp, read_tmp, nucleotides = left_nt+right_nt;
    stringstream ss;

    vector<string> nt_vec;
    list<class_snp> active_snps;
    list<class_snp>::iterator asnps_it;
    set<string>::iterator  set_it;
    map<int,int>::iterator map_it;
    map<string, class_gtf>::iterator exon_it;

    str_length = nucleotides.size();

    // Keep nearest exons to leftmost/rightmost
    if( !left_exons_map.empty() )
      {
      left_exon = (*left_exons_map.begin()).second;
      for( exon_it=left_exons_map.begin(); exon_it!=left_exons_map.end(); ++exon_it )
        if( (*exon_it).second.fin > left_exon.fin ) left_exon = (*exon_it).second;
      exon_limit=1;
      }
    if( !right_exons_map.empty() )
      {
      right_exon = (*right_exons_map.begin()).second;
      for( exon_it=right_exons_map.begin(); exon_it!=right_exons_map.end(); ++exon_it )
        if( (*exon_it).second.ini < right_exon.ini ) right_exon = (*exon_it).second;
      if(exon_limit == 0) exon_limit=2;
      else exon_limit=3;
      }
    if( exon_limit==3 && right_exon.length+left_exon.length < read_length ) exon_limit=4;

    for(unsigned k=0; k<(str_length-read_length); k++)
      {
      for( asnps_it = active_snps.begin(); asnps_it != active_snps.end(); )  // Remove past SNPs
        { if( (*asnps_it).pose <= ini_pose+k ) asnps_it = active_snps.erase( asnps_it );
          else ++asnps_it; }

      while( vec_indx < vec_len &&
             ( snp_vec[vec_indx].pose <= k+ini_pose+read_length ||  // SNPs overlapping leftmost
               snp_vec[vec_indx].pose <= k+rightmost ) )  // overlapping rightmost
        { active_snps.push_back(snp_vec[vec_indx]); vec_indx++; }
      asnps_size = active_snps.size();

      nt_vec.clear();
      nt_vec.push_back( nucleotides.substr(k,read_length) );  // This is the reference read

      if( asnps_size <= read_length*acceptable_snp_ratio && asnps_size != 0 ) // && contained_in_exon )
        {  // Process reads if less than XX% of the read is made of SNPs
        for( asnps_it = active_snps.begin(); asnps_it != active_snps.end(); ++asnps_it )
          {  // Create all alternative reads, looping through the active snps in the fragment
          nt_vec_current_size = nt_vec.size();
          for(unsigned j=0; j<nt_vec_current_size; j++)
            {
            nt_tmp = nt_vec[j];
            if( (*asnps_it).pose <= leftmost )
              nt_tmp[ (*asnps_it).pose-ini_pose-1-k ] = (*asnps_it).alt;  // leftmost SNP
            else nt_tmp[ (*asnps_it).pose-rightmost+read_length-1-k ] = (*asnps_it).alt;  // rightmost SNP
            nt_vec.push_back( nt_tmp );
            }
          }

        nt_vec_current_size = nt_vec.size();

        for(unsigned j=0; j<nt_vec_current_size; j++)
          if( exon_limit==0 ||
              (exon_limit==1 && k>=read_length-left_exon.length ) ||  // left is short
              (exon_limit==2 && k<right_exon.length-1) ||             // right is short
              (exon_limit==3 && k>=read_length-left_exon.length && k<right_exon.length-1) )
            {  // Write reads to output in fasta phred33 format
            ss.str(""); read_tmp = "";
            // Header del read, conte la info d'on s'ha generat, el num d'SNPs solapants i num del read
            ss << "@" << chr << ":" << ini_pose+k+1 << "-" << rightmost+k+1 << "j:"
              << pow(2,asnps_size) << ":" << j;
            read_tmp = ss.str() + "\n" + nt_vec[j] + "\n+\n";
            for(unsigned i=0; i<read_length; i++) read_tmp += "I";  // Quality, almost max in phred33

            out.file << read_tmp << "\n"; read_counter++;
            }
        }
      else if (asnps_size == 0 && ( exon_limit==0 ||
                (exon_limit==1 && k>=read_length-left_exon.length ) ||
                (exon_limit==2 && k<right_exon.length-1) ||
                (exon_limit==3 && k>=read_length-left_exon.length && k<right_exon.length-1) ) )
        {
        // Write reads to output in fasta phred33 format
        ss.str(""); read_tmp = "";
        // Header del read, conte la info d'on s'ha generat, el num d'SNPs solapants i num del read
        ss << "@" << chr << ":" << ini_pose+k+1 << "-" << rightmost+k+1 << "j:"
          << pow(2,asnps_size) << ":0";
        read_tmp = ss.str() + "\n" + nt_vec[0] + "\n+\n";
        for(unsigned i=0; i<read_length; i++) read_tmp += "I";  // Quality, almost max in phred33

        out.file << read_tmp << "\n"; read_counter++;
        }
      }
    }
};

/*
 *   SECONDARY FUNCTIONS
 */

static bool compare_junctp_by_leftmost( const class_junction* jp1, const class_junction* jp2 )
{ return (*jp1).leftmost < (*jp2).leftmost; }

static bool compare_junctp_by_rightmost( const class_junction* jp1, const class_junction* jp2 )
{ return (*jp1).rightmost < (*jp2).rightmost; }

bool compara_regions( class_gtf& r1, class_gtf& r2 )
{ return r1.ini < r2.ini; }

string funct_to_uppercase( string str )
{  // Convert 'acgtn' to UPPER
stringstream ss;
for(unsigned k=0; k<str.length(); k++ )
  switch( str[k] )
    {
    case 'a': ss << 'A'; break;
    case 'c': ss << 'C'; break;
    case 't': ss << 'T'; break;
    case 'g': ss << 'G'; break;
    case 'n': ss << 'N'; break;
    default: ss << str[k]; break;
    }

return ss.str();
}

void funct_write_histogram( map<int, int>& snp_density_map, class_fstream& hist_out )
{  // Outputs the histogram
map<int, int>::iterator map_it;

for( map_it = snp_density_map.begin(); map_it != snp_density_map.end(); ++map_it )
  hist_out.file << (*map_it).first << "\t" << (*map_it).second << "\n";
}


/*   =====================
 *   === MAIN FUNCTION ===
 *   =====================
 */

int main( int argc, char* argv[] )
{
unsigned chrn=25, line_count=0, vecsize=0, vecindx=0;
char line[65535];
string sys, sline, ginfopath, projpath, fastpath, snppath, genopath,
  outpath, chrt, cromosomes[24] = {
  "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
  "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
  "chr18","chr19", "chr20","chr21","chr22","chrX","chrY"};
stringstream ss;

class_gtf gtf, left_gtf, right_gtf;
class_junction junct;
class_snp snp;
class_newread newread;

list<class_newread> newread_list;
list<class_junction*> active_left_junctp_list, active_right_junctp_list;
list<class_gtf>::iterator slow_bed_it, fast_bed_it, slow_bed_it2, fast_bed_it2;
vector< list<class_gtf> > gtf_superlist(22);
vector< vector<class_snp> > snp_supervec(22);
vector< vector<class_junction> > junction_supervec(22);
vector< list<class_junction*> > junctp_byleft_superlist(22), junctp_byright_superlist(22);
list<class_junction*>::iterator left_jp_it, right_jp_it;
map<int, int> snp_density_map;
map<string, class_junction> junct_map;

class_fstream snps_in, hg19_in, gtf_in, fastq_out, hist_out, forbidden_out, count_out;

if(!debug) cout << "Starting ARTREAD_MAKER\n";
else cout << "Starting ARTREAD_MAKER -- in DEBUG mode --\n";

{  // Open all inputs and outputs
cout << "Opening all inputs and outputs" << endl;
genopath = "/home/ignasi/Imperial/Genomic_Info/";

if(!debug && argc == 5)
  {
  snppath   = argv[1];
  ginfopath = argv[2];
  outpath   = argv[3];
  read_length  = atoi(argv[4]);
  max_num_snps = atoi(argv[5]);
  }
else if(!debug && argc > 5) {cout << "Wrong number of arguments. Terminating!"; while( !cin.get() ); return 1;}
else if(debug)
  {
  ginfopath = "/home/ender/Dades/Genomic_Info/";
  projpath  = "/home/ender/projects/t2dsys/";

  snppath = projpath + "genotypes/maf/";
  outpath = projpath + "case/artificial/";

  read_length = 100;
  max_num_snps = 5;
  }
acceptable_snp_ratio = max_num_snps/(double)read_length+0.00001;
ss << read_length;

cout << "Results will be located in " << outpath << "\n";

if( !hg19_in.open   ( ginfopath+"hg19/hg19.fa", "in", "hg19" ) &&
    !gtf_in.open    ( ginfopath+"Genes_UCSC_RefSeq_Ensembl_lncRNAsNNak_hg19_noblank.gtf", "in", "gtf_file" ) &&  // With sed 's/ /_/g'

    !hist_out.open  ( outpath +"SNP_density_histogram_"+ss.str()+".txt", "out" ) &&
    !fastq_out.open ( outpath +"Artificial_reads_"+ss.str()+".fq", "out" ) &&
    !forbidden_out.open( outpath+"Polymorphic_regions_"+ss.str()+".bed", "out" ) &&
    !count_out.open  ( outpath+"generated_reads_"+ss.str()+".txt", "out" ) )
  { ; }
else { return 1; }
}

cout << "Extracting all junctions from GTF file and loading to RAM\n";
chrt=""; gtf.get_gtf( gtf_in );
while( !gtf_in.file.eof() )
  {  // Loads all gtf regions to the vector.
  if( gtf.chr != chrt )
    {
    chrt = gtf.chr;
    for(chrn=0; chrn<22; chrn++) if(chrt == cromosomes[chrn]) break;
    }

  if( chrn < 22 && gtf.type == "exon" )  // Only push exons
    gtf_superlist[chrn].push_back( gtf );

  gtf.get_gtf( gtf_in );
  }

for(chrn=0; chrn<22; chrn++)
  {  // Gets all unique junctions, and sorts two pointer lists by leftmost or rightmost
  list<class_gtf>::iterator left_it, right_it;
  left_it = gtf_superlist[chrn].begin();
  right_it = left_it; ++right_it;

  for( ; right_it != gtf_superlist[chrn].end(); ++right_it )
    {
    if( (*left_it).name == (*right_it).name )  // Both are already 'exons'
      {
      junct.make_junction( *left_it, *right_it );
      ss.str(""); ss << junct.chr << ":" << junct.leftmost << "-" << junct.rightmost;
      junct_map[ ss.str() ] = junct;  // Eliminates repeated junctions
      }

    left_it = right_it;
    }

  for( map<string,class_junction>::iterator map_it = junct_map.begin();
       map_it != junct_map.end(); ++map_it )
    junction_supervec[chrn].push_back( (*map_it).second );
  junct_map.clear();

  // Needs a separate loop since push_back reallocates memory and pointers get FUCKED UP!
  int vecsz = junction_supervec[chrn].size();
  for(int k=0; k<vecsz; k++)
    junctp_byleft_superlist[chrn].push_back ( &junction_supervec[chrn][k] );
  junctp_byright_superlist[chrn] = junctp_byleft_superlist[chrn];

  junctp_byleft_superlist[chrn].sort ( compare_junctp_by_leftmost );
  junctp_byright_superlist[chrn].sort( compare_junctp_by_rightmost );  // Finally
  }

cout << "Loading all SNPs and exons to RAM\n";

chrt="";
for( chrn=0; chrn<22; chrn++ )
  {  // Iterate over all MAF files, select binom snps at good MAF
  if( !snps_in.open( snppath+"MAF_"+cromosomes[chrn]+".frq", "in", cromosomes[chrn] ) )
    {
    snps_in.file.getline(line, 65535);  // skip the header

    snp.get_snp( snps_in );
    while( !snps_in.file.eof() )
      {  // Loads all SNPs to the RAM
      if( snp.chr != chrt )
        {
        chrt = snp.chr;
        for(chrn=0; chrn<22; chrn++) if(chrt == cromosomes[chrn]) break;
        }

      if(chrn < 22 && snp.isgood() )
        {
        snp.make_nt_chars();
        snp_supervec[chrn].push_back(snp);
        }
      snp.get_snp( snps_in );
      }

    snps_in.close();
    }
  else return 1;

  cout << "Read " << snp_supervec[chrn].size() << " SNPs for "
    << cromosomes[chrn] << "\n";
  }


cout << "Starting read creation\n";

hg19_in.file.getline(line, 65535); sline = line;
while( !hg19_in.file.eof() )
  {
  if(line[0] == '>')
    {  // Change of chr
    string strline(line);
    chrt = strline.substr(1, strline.find('\0') );
    for(chrn=0; chrn<22; chrn++) if(cromosomes[chrn] == chrt) break;
    cout << "Working through " << chrt << "\nWe have written "
      << read_counter << " reads so far\n";

    if(chrn < 22)
      {  // Init
      vecsize = snp_supervec[chrn].size();
      gtf_superlist[chrn].sort( compara_regions );
      slow_bed_it = gtf_superlist[chrn].begin();
      slow_bed_it2 = gtf_superlist[chrn].begin();
      line_count = 0; vecindx = 0; snp = snp_supervec[chrn][vecindx];
      left_jp_it  = junctp_byleft_superlist[chrn].begin();
      right_jp_it = junctp_byright_superlist[chrn].begin();
      }
    }

  if(chrn<22)  // Since we are ignoring chrX and chrY for this project
    {
    // Adds the current nucleotides to all newreads in tghe list, and writes the finalized ones
    if( !newread_list.empty() )
      for( list<class_newread>::iterator it=newread_list.begin(); it!=newread_list.end(); )
        {
        // Add to newread all exons that start in this line
        for( fast_bed_it = slow_bed_it;
             ( fast_bed_it != gtf_superlist[chrn].end() ) &&
             ( (*fast_bed_it).ini <= line_count*line_char_count ) &&
             ( (*fast_bed_it).fin >= (line_count-1)*line_char_count+1 ) ; ++fast_bed_it )
          (*it).exons_map[ (*fast_bed_it).key ] = *fast_bed_it;

        (*it).update( sline );  // Updates the nt string

        if( (*it).status == 1 )
          {  // Creates all reads and writes them out, along with the too-polymorphic regions
          (*it).write_reads( fastq_out, forbidden_out, snp_density_map );
          it = newread_list.erase(it);  // Already advances the iterator
          }
        else ++it;
        }
    while( slow_bed_it != gtf_superlist[chrn].end() &&
           line_count*line_char_count > (*slow_bed_it).fin  )
      ++slow_bed_it;  // Slow iteration, to avoid skipping overlapping or long exons

    // Iniciating, updating and writting junction reads
    while( left_jp_it != junctp_byleft_superlist[chrn].end() &&
           (*left_jp_it)->leftmost-read_length+1 <= line_count*line_char_count )
      {  // Initiates the leftmost part of the junction if NO SNP overlaps
      active_left_junctp_list.push_back( *left_jp_it );
      ++left_jp_it;  // Advances the left jp interator if SNP does not overlap
      }
    while( right_jp_it != junctp_byright_superlist[chrn].end() &&
           (*right_jp_it)->rightmost <= line_count*line_char_count )
      {  // Initiates the rightmost part of the junction if NO SNP overlaps
      active_right_junctp_list.push_back( *right_jp_it );
      ++right_jp_it;  // Advances the right jp interator if SNP does not overlap
      }

    if( !active_left_junctp_list.empty() )
      for( list<class_junction*>::iterator ajp_it = active_left_junctp_list.begin();
           ajp_it!=active_left_junctp_list.end(); ++ajp_it )
        {  // Iterates the list of active left junctions, and updates or closes them
        // Add to left junction all exons starting in this line
        for( fast_bed_it2 = slow_bed_it2;
             ( fast_bed_it2 != gtf_superlist[chrn].end() ) &&
             ( (*fast_bed_it2).ini <= line_count*line_char_count ) &&
             ( (*fast_bed_it2).fin >= (line_count-1)*line_char_count+1 ) ; ++fast_bed_it2 )
          if( (*fast_bed_it2).length <= read_length )
            (*ajp_it)->left_exons_map[ (*fast_bed_it2).key ] = *fast_bed_it2;

        if( (*ajp_it)->leftmost-read_length+1 <= line_count*line_char_count &&
            (*ajp_it)->leftmost >= line_count*line_char_count )
          (*ajp_it)->update_left_nt( sline, line_count );
        else if( (*ajp_it)->leftmost-read_length+1 <= line_count*line_char_count &&
                 (*ajp_it)->leftmost < line_count*line_char_count )
          (*ajp_it)->close_left = 1;
        }
    if( !active_right_junctp_list.empty() )
      for( list<class_junction*>::iterator ajp_it = active_right_junctp_list.begin();
           ajp_it!=active_right_junctp_list.end(); ++ajp_it )
        {  // Iterates the list of active right junctions, and updates or closes them
        // Add to right junction all exons starting in this line
        for( fast_bed_it2 = slow_bed_it2;
             ( fast_bed_it2 != gtf_superlist[chrn].end() ) &&
             ( (*fast_bed_it2).ini <= line_count*line_char_count ) &&
             ( (*fast_bed_it2).fin >= (line_count-1)*line_char_count+1 ) ; ++fast_bed_it2 )
          if( (*fast_bed_it2).length <= read_length )
            (*ajp_it)->right_exons_map[ (*fast_bed_it2).key ] = *fast_bed_it2;

        if( (*ajp_it)->rightmost <= line_count*line_char_count &&
            (*ajp_it)->rightmost+read_length-1 >= line_count*line_char_count )
          (*ajp_it)->update_right_nt( sline, line_count );
        else if( (*ajp_it)->rightmost < line_count*line_char_count &&
                 (*ajp_it)->rightmost+read_length-1 < line_count*line_char_count )
          (*ajp_it)->close_right=1;
        }
    while( slow_bed_it2 != gtf_superlist[chrn].end() &&
           line_count*line_char_count > (*slow_bed_it2).fin  )
      ++slow_bed_it2;  // Slow iteration, to avoid skipping overlapping or long exons

    // Adds the overlapping SNPs to newreads or junctreads
    while( vecindx < vecsize && line_count*line_char_count > snp.pose-read_length )
      {  // Init newreads and update the SNPs of those already init
      if( !newread_list.empty() &&
          snp.pose - newread_list.back().snp_vec.back().pose <= read_length )
        newread_list.back().snp_vec.push_back(snp);  // Add the SNP if <read_length from previous
      else
        {  // If list is empty or next SNP is too far, make a new one
        newread.create( snp, line_count, sline );
        newread_list.push_back(newread);
        }

      while( left_jp_it != junctp_byleft_superlist[chrn].end() &&
             (*left_jp_it)->leftmost-read_length+2 <= snp.pose )
        {  // Initiates the leftmost part of the junction if a SNP overlaps
        active_left_junctp_list.push_back( *left_jp_it );
        ++left_jp_it;  // Advances the left jp interator if SNP overlaps
        }
      while( right_jp_it != junctp_byright_superlist[chrn].end() &&
             (*right_jp_it)->rightmost <= snp.pose )
        {  // Initiates the rightmost part of the junction if a SNP overlaps
        active_right_junctp_list.push_back( *right_jp_it );
        ++right_jp_it;  // Advances the right jp interator if SNP overlaps
        }

      if( !active_left_junctp_list.empty() )
        for( list<class_junction*>::iterator ajp_it = active_left_junctp_list.begin();
             ajp_it!=active_left_junctp_list.end(); ++ajp_it )
          {  // Iterates the list of active left junctions, and adds overlapping snps
          if( (*ajp_it)->leftmost-read_length+2 <= snp.pose &&
              (*ajp_it)->leftmost >= snp.pose )
            (*ajp_it)->snp_vec.push_back( snp );
          }
      if( !active_right_junctp_list.empty() )
        for( list<class_junction*>::iterator ajp_it = active_right_junctp_list.begin();
             ajp_it!=active_right_junctp_list.end(); ++ajp_it )
          {  // Iterates the list of active right junctions, and adds overlapping snps
          if( (*ajp_it)->rightmost <= snp.pose &&
              (*ajp_it)->rightmost+read_length-2 >= snp.pose )
            (*ajp_it)->snp_vec.push_back( snp );
          }

      vecindx++; // Get next SNP
      if(vecindx < vecsize) snp = snp_supervec[chrn][vecindx];
      else break;
      }

    // Performing the actual junction read closure
    if( !active_left_junctp_list.empty() )
      for( list<class_junction*>::iterator ajp_it = active_left_junctp_list.begin();
           ajp_it!=active_left_junctp_list.end(); )
        { if( (*ajp_it)->close_left )
            { (*ajp_it)->close_left_nt( sline );
              ajp_it = active_left_junctp_list.erase(ajp_it); }
          else ++ajp_it; }

    if( !active_right_junctp_list.empty() )
      for( list<class_junction*>::iterator ajp_it = active_right_junctp_list.begin();
           ajp_it!=active_right_junctp_list.end(); )
        { if( (*ajp_it)->close_right )
            { (*ajp_it)->close_right_nt( sline );
              (*ajp_it)->write( fastq_out );  // Write the junctions when rightmost is closed
              ajp_it = active_right_junctp_list.erase(ajp_it); }
          else ++ajp_it; }
    }

  hg19_in.file.getline(line, 65535); sline = line;
  line_count++;
  }

// Write the densty of SNPs histogram
funct_write_histogram( snp_density_map, hist_out );

cout << "\nWe have written a total of " << read_counter << " reads\n";
count_out.file << read_counter << endl;

cout << "\nAll done!\n";

snps_in.close(); hg19_in.close(); fastq_out.close(); gtf_in.close();
hist_out.close(); forbidden_out.close(); count_out.close();

return 0;
}














