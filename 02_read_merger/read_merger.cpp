/*
 *  READ_MERGER.CPP
 *
 *  Reads an accepted_hits.sorted.sam ENHANCED and MASKED alignments,
 *  and outputs a MERGED accepted_hits.sorted.bam and non-clonal
 *  accepted_hits.sorted.bam files.
 *
 *  Written by Ignasi Moran
 */

#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <map>
#include <set>
#include <iomanip>
#include "../utils/class_fstream.h"


using namespace std;
const bool debug = 0,  // DEBUG FLAG
  account_for_discordant = 1,  // Account for possible chr-discordant alignments
  write_intermediate = 0;  // Writes separated common, enhSp and maskSp sam files

const string chrequal = "=", chrdummy = "chrDummy", cromosomes[24] = {
  "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
  "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
  "chr18","chr19", "chr20","chr21","chr22","chrX","chrY"};

/*
 *   CLASSES
 */

class class_read
{  // All fields related to an RNAseq read
public:
  bool found;
  int flag, score, align_length;
  long int pose;
  string name, cigar, quality, mate_pose, additional,
    key;  // search key = name{bases , uses { as separator since it is not found in quality characters
  const string *chr, *mate_chr;

  class_read() : found(false), flag(0), score(0), align_length(0), pose(0),
     name(""), cigar(""), quality(""), mate_pose(""), additional("."),
     chr( &chrdummy ), mate_chr( &chrdummy ) {}  // Default constructor to back our asses

  void get( char l[1024], bool noncl=0 )
    {  // Parses all the sam info into variables from the whole getline, and generates a key
    int num_sam_camps = 12;
    size_t str_pose, prev_pose;
    string dummy, line(l);  // string conversion
    stringstream noncl_key_ss;  // Which is "pose{cigar{mate_pose{bases"

    prev_pose = 0;
    str_pose = line.find("\t", prev_pose);
    for(int k=1; k<num_sam_camps && str_pose != string::npos; k++)
      {
      dummy = line.substr(prev_pose, str_pose-prev_pose);

      switch(k)
        {  // Reads the read info and stores it to write it out later
        case 1: if( !noncl ) key = dummy;
          else name = dummy; break;
        case 2: flag = atoi(dummy.c_str()); break;
        case 3:  // To reduce ram footprint, chrs are stored as pointers
          chr = &chrdummy;
          for( int chrn=0; chrn<24; chrn++ )
            if( cromosomes[chrn] == dummy ) chr = &(cromosomes[chrn]);
          break;
        case 4: if(!noncl) pose = atol(dummy.c_str());
          else noncl_key_ss << dummy << "{"; break;
        case 5: score = atoi(dummy.c_str()); break;
        case 6: if(!noncl) cigar = dummy;
          else noncl_key_ss << dummy << "{"; break;
        case 7:
          mate_chr = &chrdummy;
          if( dummy == "=" ) mate_chr = &(chrequal);
          else{ for( int chrn=0; chrn<24; chrn++ )
            if( cromosomes[chrn] == dummy ) mate_chr = &(cromosomes[chrn]); }
          break;
        case 8: if(!noncl) mate_pose = dummy;  // Was atol, but it is not used as numeric
          else noncl_key_ss << dummy << "{"; break;
        case 9: align_length = atoi(dummy.c_str()); break;
        case 10:
          if( !noncl ) key = key + "{" + dummy;
          else noncl_key_ss << dummy; break;
        case 11: quality = dummy; break;
        default: break;
        }

      prev_pose = str_pose+1;
      str_pose = line.find("\t", prev_pose);
      }

    // Complex key to remove clonal reads: if they are equal in all
    // these fields they will be discarded
    if( noncl ) key = noncl_key_ss.str();

    // Additional fields, stored all in this variable
//    str_pose = line.find("XR:Z:"); str_pose--;  // Temporary fix to samtools parse error
//    if( str_pose != string::npos ) additional = line.substr(prev_pose, str_pose-prev_pose);  // If enhanced
//    else
    additional = line.substr(prev_pose, line.size()-prev_pose);  //If masked
    }

  string print( bool noncl=0 )
    {  // Returns a string to print read in sam format
    stringstream ss;
    size_t p, p2;
    string ps, cg, matep, bp;

    if( noncl == 0 )
      {  // key is "name{bases"
      ss << key.substr( 0, key.find("{") ) << "\t"
        << flag << "\t" << *chr << "\t" << pose << "\t" << score << "\t"
        << cigar << "\t" << *mate_chr << "\t" << mate_pose << "\t" << align_length << "\t"
        << key.substr( key.find("{")+1, string::npos ) << "\t"
        << quality << "\t" << additional;
      }

    else
      {  // key is "pose{cigar{mate_pose{bases"
      p=0; for(int i=0; i<3; i++)
        {
        p = key.find("{", p)+1;
        if( i==0 ) { ps = key.substr(0, p-1); p2=p; }
        if( i==1 ) { cg = key.substr( p2, p-p2-1 ); p2=p; }
        if( i==2 ) { matep = key.substr( p2, p-p2-1 ); bp = key.substr( p, string::npos ); }
        }

      ss << name << "\t" << flag << "\t" << *chr << "\t" << ps << "\t" << score << "\t"
        << cg << "\t" << *mate_chr << "\t" << matep << "\t" << align_length << "\t"
        << bp << "\t" << quality << "\t" << additional;
      }

    return ss.str();
    }
};

/*
 *   SECONDARY FUNCTIONS
 */

int funct_search_read_in_map( class_read& enh_read, map<string,class_read>& mask_reads_map,
  map<string,class_read>::iterator& map_it )
{  // Searches for the read in the map, and returns 1, 2 or 3 depending if found and if specific
class_read dummy_read;

map_it = mask_reads_map.find( enh_read.key );  // Searches the enh read in the masked map

if( map_it != mask_reads_map.end() )
  {  // Fount in the map
  (*map_it).second.found = true;  // Marked as found, to discard as mask-sp

  if( (*map_it).second.cigar.find("N") == string::npos )
    return 1;  // Report a match without splicing
  else return 2;  // Report a match with splicing
  }
else return 3;  // Report no match in the masked map

return 0;  // Defaults
}

int funct_match_read_ends( class_read& mread, class_read& eread )
{  // Returns the number of ends that are common between the two reads, either 0, 1 or 2
char delims[7] = {'M','N','S','I','D','H','P'};  // M normal align, N intron, S soft clip, I insertion, D deletion
size_t sposem, mposem, sposee, mposee, str_pose1, str_pose2, min_pose;
int num_matches=0;
long int fixed_inim=mread.pose, fixed_inie=eread.pose, calc_finm, calc_fine;

sposem=mread.cigar.find( delims[2] ); mposem=mread.cigar.find( delims[0] );
sposee=eread.cigar.find( delims[2] ); mposee=eread.cigar.find( delims[0] );

if( sposem != string::npos && sposem < mposem )  // S present and before M
  fixed_inim = mread.pose - atol( mread.cigar.substr( 0, sposem ).c_str() );
if( sposee != string::npos && sposee < mposee )  // S present and before M
  fixed_inie = eread.pose - atol( eread.cigar.substr( 0, sposee ).c_str() );

if( fixed_inim == fixed_inie ) num_matches=1;

// Now fun with CIGARs!! to try and calculate the rightmost end
str_pose1=0; str_pose2=0; calc_finm = fixed_inim;
while( str_pose1 != string::npos && str_pose1 < mread.cigar.length() )
  {
  min_pose = 1000;
  for( unsigned u=0; u<7; u++ )
    {  // Looks for the next non-numeric character in the CIGAR string
    str_pose2 = mread.cigar.find( delims[u], str_pose1 );
    if( str_pose2 != string::npos ) min_pose = min( min_pose, str_pose2 );
    }
  if( min_pose != 1000 )
    {  // Calculates the other end of the masked read (since enhanced will never be broken)
    calc_finm += atol( mread.cigar.substr( str_pose1, min_pose-str_pose1 ).c_str() );
    str_pose1 = ++min_pose;  // Should increment before the =, so that the while triggers.
    }
  }

// Repeat for enhanced read
str_pose1=0; str_pose2=0; calc_fine = fixed_inie;
while( str_pose1 != string::npos && str_pose1 < eread.cigar.length() )
  {
  min_pose = 1000;
  for( unsigned u=0; u<7; u++ )
    {  // Looks for the next non-numeric character in the CIGAR string
    str_pose2 = eread.cigar.find( delims[u], str_pose1 );
    if( str_pose2 != string::npos ) min_pose = min( min_pose, str_pose2 );
    }
  if( min_pose != 1000 )
    {  // Calculates the other end of the masked read (since enhanced will never be broken)
    calc_fine += atol( eread.cigar.substr( str_pose1, min_pose-str_pose1 ).c_str() );
    str_pose1 = ++min_pose;  // Should increment before the =, so that the while triggers.
    }
  }

if( calc_finm == calc_fine ) num_matches++;

return num_matches;
}

void funct_sam2bam( string& file )
{  // Converts a given sam file to bam, sorts and indexes it
string sys, path = file.substr( 0, file.find_last_of("/")+1 );

sys = "sleep 10; samtools view -bhS " + file + " | samtools sort -@ 12 - " + path
  + "accepted_hits.sorted; samtools index " + path + "accepted_hits.sorted.bam; sleep 10;";
if(debug) cout << sys << "\n";
system(sys.c_str());
//if( system(sys.c_str()) == -1 ) { stringstream ss;
//  ss << "Sam2bam system call error with " << sys << "\n"; throw ss.str(); }

if(!debug) { sys = "rm " + file; system(sys.c_str()); }
//  if( system(sys.c_str()) == -1 ) { stringstream ss;
//    ss << "Sam2bam system call error with " << sys << "\n"; throw ss.str(); } }
}

void funct_write_stats( long int& found_count, long int& concord_count, long int& chr_discordant,
  long int& enh_count, long int& mask_count, long int& masked_read_count,
  long int& enh_read_count, long int& ns_agree, long int& ns_disagree,
  long int& s_agree, long int& s_disagree, long int& clon_count, class_fstream& stats_out )
{  // Statistics
double total = concord_count+enh_count+mask_count;
stringstream ss;

ss << setprecision(3) << fixed
  << "Statistics: Got " << masked_read_count/(double)1000000 << "M masked reads and "
  << enh_read_count/(double)1000000 << "M enhanced reads\nReported a total of "
  << total/(double)1000000 << "M unique reads (a mean of "
  << total*2*100/(double)(masked_read_count+enh_read_count)-100
  << "% increase), of which\n " << enh_count/(double)1000000 << "M enhanced specific (" << enh_count*100/total
  << "%),\n " << mask_count/(double)1000000 << "M masked specific (" << mask_count*100/total
  << "%),\n " << concord_count/(double)1000000 << "M concordant (" << concord_count*100/total
  << "%), with\n  " << ns_agree/(double)1000000 << "M nosplice-agree (" << ns_agree*100/(double)found_count
  << "%), \n  " << ns_disagree/(double)1000000 << "M nosplice-disagree (" << ns_disagree*100/(double)found_count
  << "%), \n  " << s_agree/(double)1000000 << "M splice-agree (" << s_agree*100/(double)found_count
  << "%), i\n  " << s_disagree/(double)1000000 << "M splice-disagree (" << s_disagree*100/(double)found_count
  << "%),\n " << chr_discordant/(double)1000000 << "M chromosome discordant (" << chr_discordant*100/total << "%)\n"
  << clon_count/(double)1000000 << "M non-clonal reads also writen out (" << clon_count*100/(double)total << "%)\n";

cout << "\n" << ss.str();
stats_out.file << ss.str();
}

/*   =====================
 *   === MAIN FUNCTION ===
 *   =====================
 */

int main( int argc, char* argv[] )
{
char maskline[1024], enhline[1024], headline[1024];
long int masked_read_count = 0, merged_read_count = 0, s_agree=0, s_disagree=0, ns_agree=0, ns_disagree=0,
  found_count=0, concord_count=0, enh_read_count=0, enh_count=0, mask_count=0, chr_discordant=0, clon_count=0;  // Counters for general statistics

string sys, header_file, enhanced_file, masked_file,
  out_standard_path, out_nonclonal_path, chr_t ="";
stringstream header_ss;

vector<class_read> enhancedsp_vector;
map<string, class_read> mask_reads_map, specific_reads_map;
map<string, class_read>::iterator map_it;

class_fstream config_in, header_in, enhanced_in, masked_in, merged_out, nonclon_out, stats_out,
  enhanced_out, masked_out, merged_in, common_out, discord_out;

if(!debug) cout << "Initiating READ MERGER\n";
else cout << "Initiating READ MERGER -- Debug mode --\n";

{  // Opening inputs and outputs
unsigned num_arguments = 4;
if(!debug)
  {  // Read input arguments from config file
  if( !config_in.open( "config.ini", "in", "config" ) )
    {
    string config_argv[num_arguments+1], dummy;
    cout << "Read from config file:\n";
    for(unsigned k=1; k<num_arguments+1; k++){
      config_in.file >> dummy >> config_argv[k];
      cout << dummy <<"\t"<< config_argv[k] << "\n";}

    enhanced_file = config_argv[1];
    masked_file = config_argv[2];
    out_standard_path = config_argv[3];
    out_nonclonal_path = config_argv[4];
    }
  else { cout << "Error reading the config file, terminating!"; return 1; }
  config_in.close();
  }
else if(debug && argc == (int)num_arguments+1)
  {  // In debug mode, use command line arguments instead
  enhanced_file = argv[1];
  masked_file = argv[2];
  out_standard_path = argv[3];
  out_nonclonal_path = argv[4];
  }
else if(debug && argc != (int)num_arguments+1)
  {cout << "Wrong number of arguments. Terminating!"; while( !cin.get() ); return 1;}

cout << "Given paths were:\nEnhanced: " << enhanced_file
  << "\nMasked: " << masked_file
  << "\nOutputs: " << out_standard_path << " and " << out_nonclonal_path << "\n";

sys = "mkdir " + out_standard_path + "; mkdir " + out_nonclonal_path;
system(sys.c_str());
//if( system(sys.c_str()) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}

// Extracts the header for further use
header_file = out_standard_path+"header.sam";
sys = "sleep 10; samtools view -H " + masked_file + " > "
  + header_file + "; sleep 10;";
system(sys.c_str());
//if( system(sys.c_str()) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}

if( !header_in.open   ( header_file, "in", "header") &&
    !merged_out.open  ( out_standard_path +"accepted_hits.merged.sam", "out" ) &&
    !nonclon_out.open ( out_nonclonal_path+"accepted_hits.nonclonal.sam", "out" ) &&
    !stats_out.open   ( out_standard_path +"Merging_statistics.txt", "out" ) )
  {
  // Read the header, close it and write the output files
  int i=0; header_in.file.getline(headline, 1024);
  while( !header_in.file.eof() )
    {
    if(i<25) header_ss << headline << "\n";
    header_in.file.getline(headline, 1024);
    i++;
    }
  header_in.close(); sys = "rm " + header_file;
  system(sys.c_str());
//  if( system(sys.c_str()) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}

  // Writing output headers
  merged_out.file << header_ss.str() << "@PG\tID:Bowtie2_and_TopHat2_merged\n";
  nonclon_out.file << header_ss.str() << "@PG\tID:Bowtie2_and_TopHat2_merged_nonclonal\n";
  }
else { while( !cin.get() ) ; return 1; }

if( write_intermediate )
  {  // Intermediate files
  if( !enhanced_out.open( out_standard_path+"accepted_hits.enhanced_specific.sam", "out" ) &&
      !masked_out.open  ( out_standard_path+"accepted_hits.masked_specific.sam", "out" ) &&
      !common_out.open  ( out_standard_path+"accepted_hits.common.sam", "out" ) &&
      !discord_out.open ( out_standard_path+"accepted_hits.discordant.sam", "out" ) )
    {
    enhanced_out.file << header_ss.str() << "@PG\tID:Bowtie2_enhanced_specific\n";
    masked_out.file << header_ss.str()   << "@PG\tID:TopHat2_masked_specific\n";
    common_out.file << header_ss.str()   << "@PG\tID:Bowtie2_and_TopHat2_common\n";
    discord_out.file << header_ss.str()  << "@PG\tID:Bowtie2_and_TopHat2_common\n";
    }
  }
}

string mask_temp_file = out_standard_path+"mask_temp.sam",
  enh_temp_file = out_standard_path+"mask_temp.sam";
class_read enh_read, mask_read;

//int enhcounter=0;
for(int chrn=0; chrn<24; chrn++) if( !debug || chrn == 19 || chrn == 20 || chrn == 21 )
  {  // Chromosome loop
  chr_t = cromosomes[chrn];
  cout << "\nStart of " << chr_t << ". Extracting the reads from the bam files\n";

  // Extract the bam files reads, since they are randomly sorted
  sys = "sleep 10; samtools view " + masked_file + " " + chr_t + " > "
    + mask_temp_file + "; sleep 10;";
  system(sys.c_str());  // Unprotected calls, to see if we get more info when they blow up
//  if( system(sys.c_str())  == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}
  sys = "sleep 10; samtools view " + enhanced_file + " " + chr_t + " > "
    + enh_temp_file + "; sleep 10;";
  system(sys.c_str());  // Unprotected calls, to see if we get more info when they blow up
//  if( system(sys.c_str())  == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}

  if( !enhanced_in.open ( enh_temp_file, "in", chr_t+" enhanced reads" ) &&
      !masked_in.open   ( mask_temp_file, "in", chr_t+" masked reads" ) )
    { enhanced_in.file.getline(enhline, 1024); enh_read.get( enhline );
      masked_in.file.getline(maskline, 1024); mask_read.get( maskline ); }
    else { while( !cin.get() ) ; return 1; }

  cout << "Creating a map of the masked reads\n";
  while( !masked_in.file.eof() && chr_t == *(mask_read.chr) )
    {  // Reads the masked reads and makes a map
    mask_read.get( maskline );
    mask_reads_map[ mask_read.key ] = mask_read;  // This destroys some name+basepair duplicates!

    masked_read_count++; if( masked_read_count%1000000 == 0 )
      cout << "Processed " << masked_read_count/1000000 << "M masked reads so far\n";

    masked_in.file.getline(maskline, 1024);
    }

  cout << "Comparing the enhanced reads against the masked map\n";

  while( !enhanced_in.file.eof() && chr_t == *(enh_read.chr) )  // Assuming same ordering for mask and enh
    {  // Reads the enhanced reads and compares it to the map
    enh_read.get( enhline );
//    enhcounter++; cout << enhcounter << " " << enhline << "\n"; cin.get();
    switch( funct_search_read_in_map( enh_read, mask_reads_map, map_it ) )
      {  // Searches the enh read in the masked map, deciding the outcome depending on the match found
      case 1:  // Found a match in the masked map and it was not spliced
        if( write_intermediate ) common_out.file << enh_read.print() << "\n";
        if( enh_read.chr == (*map_it).second.chr &&
            funct_match_read_ends( (*map_it).second, enh_read ) == 2 )  // The 2 ends match
          { ns_agree++; merged_out.file << enh_read.print() << "\n"; }  // Prints bowtie read
        else { ns_disagree++; if( write_intermediate )  // Also prints both discordant reads
          discord_out.file << enh_read.print() << "\n" << (*map_it).second.print() << "\n"; }
        break;
      case 2:  // Found a match that was spliced - check for one end consistency between bowtie and star
        if( write_intermediate ) common_out.file << (*map_it).second.print() << "\n";
        if( enh_read.chr == (*map_it).second.chr &&
            funct_match_read_ends( (*map_it).second, enh_read ) > 0 )  // At least one end matches
          { s_agree++; merged_out.file << (*map_it).second.print() << "\n"; }  // Prints star read
        else { s_disagree++; if( write_intermediate )  // Also prints both discordant reads
          discord_out.file << enh_read.print() << "\n" << (*map_it).second.print() << "\n"; }
        break;
      case 3:  // Match not found: likely enhanced specific read
        if( account_for_discordant )
          enhancedsp_vector.push_back( enh_read );  // For future comparison with rest of chr
        else
          {  // If low RAM: Assume it is enh specific
          if( write_intermediate ) enhanced_out.file << enh_read.print() << "\n";
          merged_out.file << enh_read.print() << "\n";
          enh_count++;
          }
        break;
      default: cout << "Search read switch just defaulted!\n"; cin.get(); break;
      }

    enh_read_count++;
    found_count = s_agree+s_disagree+ns_agree+ns_disagree;
    concord_count = s_agree+ns_agree;
    if( enh_read_count%1000000 == 0 )
      cout << "Compared " << enh_read_count/1000000
        << "M enhanced reads so far\n";

    enhanced_in.file.getline(enhline, 1024);
    }

  for(map_it = mask_reads_map.begin(); map_it != mask_reads_map.end(); ++map_it)
    if( !(*map_it).second.found )
      {
      if( account_for_discordant )  // Searches for masked specific reads, and writes them in their own map
        specific_reads_map[ (*map_it).second.key ] = (*map_it).second;  // Fills multi-chr masksp map
      else  // Not accounting for chr-discordant reads
        {
        if( write_intermediate ) masked_out.file << (*map_it).second.print() << "\n";
        merged_out.file << (*map_it).second.print() << "\n";
        mask_count++;
        }
      }
  mask_reads_map.clear();

  cout << "End of " << chr_t << "\n";
  enhanced_in.close(); masked_in.close();
  }

if( account_for_discordant )
  {  // End moves: Dealing with chr-discordant 'specific' reads
  int enh_size = enhancedsp_vector.size();
  cout << "\nResolving chr-discordant potentially specific reads\n"
    << "Enhanced specific vector contains " << enh_size << " elements,\n"
    << "masked specific map contains " << specific_reads_map.size() << "\n";

  for(int j=0; j<enh_size; j++)
    {
    switch( funct_search_read_in_map( enhancedsp_vector[j], specific_reads_map, map_it ) )
      {  // Searches the 'masked specific' map against the 'enhanced specific' reads
      case 1:  // Found a non spliced match in the masked map, but is in a diff chromosome
      case 2:  // Found a match that was spliced - the read was aligned to diff chrs, so we discard it
        if( write_intermediate ) common_out.file << (*map_it).second.print() << "\n";
        chr_discordant++;
        break;
      case 3:  // Match not found: enhanced specific read
        if( write_intermediate ) enhanced_out.file << enhancedsp_vector[j].print() << "\n";
        merged_out.file << enhancedsp_vector[j].print() << "\n";
        enh_count++;
        break;
      default: cout << "Search read switch2 just defaulted!\n"; cin.get(); break;
      }

    if( (chr_discordant+enh_count)%1000000 == 0 )
      cout << "Processed " << (chr_discordant+enh_count)/1000000
        << "M of prob enhanced-specific reads so far\n";
    }

  for(map_it = specific_reads_map.begin(); map_it != specific_reads_map.end(); ++map_it)
    if( !(*map_it).second.found )  // Remaining reads are completely masked_specific
      {
      if( write_intermediate ) masked_out.file << (*map_it).second.print() << "\n";
      merged_out.file << (*map_it).second.print() << "\n";
      mask_count++;
      }
  specific_reads_map.clear();
  }

cout << "\nConverting merged output to .bam and sorting\n";
merged_out.file.flush(); merged_out.close();
funct_sam2bam( sys = out_standard_path+"accepted_hits.merged.sam" );
//try { funct_sam2bam( sys = out_standard_path+"accepted_hits.merged.sam" ); }
//catch( string sam2bamerror ){ cout << "System error in sam2bam with "
//  << sam2bamerror << "\nTerminating!\n"; cin.get(); return 1; }

// Now we read the merged sam file, put it in a hash, and print that out, to create a nonclonal output
cout << "\nDone. Now removing clonal reads";  // Done separately rather than in parallel, to reduce ram footprint

for(int chrn=0; chrn<24; chrn++) if( !debug || chrn == 19 || chrn == 20 || chrn == 21 )
  {
  chr_t = cromosomes[chrn];
  cout << "\nStart of " << chr_t << ". Extracting the reads from the merged sam file.\n";

  // Extract the chr reads from the merged sam file
  sys = "sleep 10; samtools view " + out_standard_path + "accepted_hits.sorted.bam " + chr_t + " > "
    + out_standard_path + "temp.sam; sleep 10;";
  system(sys.c_str());
//  if( system(sys.c_str()) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}

  if( !merged_in.open ( out_standard_path +"temp.sam", "in", chr_t+" merged reads" ) )
    { merged_in.file.getline(maskline, 1024); mask_read.get( maskline,1 ); }
  else { while( !cin.get() ) ; return 1; }

  cout << "Creating a map of the merged reads\n";
  while( !merged_in.file.eof() && chr_t == *(mask_read.chr) )
    {  // Reads the masked reads and makes a map
    mask_read.get( maskline, 1 );
    mask_reads_map[ mask_read.key ] = mask_read;  // This destroys clonal duplicates

    merged_read_count++; if( merged_read_count%1000000 == 0 )
      cout << "Processed " << merged_read_count/1000000 << "M merged reads so far\n";

    merged_in.file.getline(maskline, 1024);
    }
  merged_in.close();

  cout << "Writing " << mask_reads_map.size() << " nonclonal reads\n";
  for(map_it = mask_reads_map.begin(); map_it != mask_reads_map.end(); ++map_it)
    {
    nonclon_out.file << (*map_it).second.print(1) << "\n";
    clon_count++;
    }
  mask_reads_map.clear();
  }

cout << "Converting non-clonal output to .bam and sorting\n";
nonclon_out.file.flush(); nonclon_out.close();

funct_sam2bam( sys = out_nonclonal_path+"accepted_hits.nonclonal.sam" );
//try { funct_sam2bam( sys = outnonclonal+"accepted_hits.nonclonal.sam" ); }
//catch( string sam2bamerror ){ cout << "System error in sam2bam with "
//  << sam2bamerror << "\nTerminating!\n"; cin.get(); return 1; }

// Statistics
funct_write_stats( found_count, concord_count, chr_discordant, enh_count, mask_count,
  masked_read_count, enh_read_count, ns_agree, ns_disagree, s_agree, s_disagree, clon_count, stats_out );


cout << "\nDone!\n";

sys = "rm " + enh_temp_file; system( sys.c_str() );
//if( system( sys.c_str() ) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}
sys = "rm " + mask_temp_file; system( sys.c_str() );
//if( system( sys.c_str() ) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}
sys = "rm " + out_standard_path + "temp.sam"; system( sys.c_str() );
//if( system( sys.c_str() ) == -1 ) {cout << "System error in " << sys << "\nTerminating!\n"; cin.get(); return 1;}

header_in.close(); enhanced_in.close(); masked_in.close(); nonclon_out.close();
stats_out.close(); enhanced_out.close(); masked_out.close(); common_out.close();

return 0;
}










