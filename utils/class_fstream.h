/*
 *  class_fstream.h
 *
 *  Aquest header conte class_fstream, que gestiona els arxius d'entrada i sortida awsomelyment
 *
 */

#include <iostream>
#include <fstream>
#include <string>

class class_fstream
{  // Conte el fstream i la funcio per obrirlo be (amb path = "/home/ignasi/IDIBAPS/Alignments/"+p)
public:
  std::string path, in_or_out, tag;
  std::fstream file;

  bool open(std::string p, std::string in_or_out, std::string t = "", int verbose=1)
    {  // Opens the given file path
    int pos, prepos;
    char c[255];
    std::string sys;
    FILE *pope;

    path = p;

    if( t == "") tag = p;
    else tag = t;

    if(in_or_out == "in")
      {
      if(verbose) std::cout << "Opening input " + tag + "\n";
      file.open(path.c_str(), std::ios::in);
      if( !file.is_open() )
        {
        file.close();
        std::cout << "Not possible to open " + path + ". Terminating!";
        return 1;
        }
      }
    else if(in_or_out == "out")
      {
      if(verbose) std::cout << "Opening output " + p + "\n";

      // Checking path to the output file exists
      for(pos=0; pos != -1 ; ) {prepos = pos; pos++; pos = path.find('/', pos);}
      sys = "if test -d \"" + path.substr(0, prepos) + "\"; then echo \"T\"; else echo \"F\"; fi";
      pope = popen(sys.c_str(), "r"); if(fgets(c, 255, pope) == NULL) std::cout << "Fgets ERROR!\n";

      if( c[0] == 'T') file.open(path.c_str(), std::ios::out);
      else if( c[0] == 'F')
        {
        std::cout << "Non existing path to " + path + ". Terminating!";
        return 1;
        }
      else {std::cout << "class_fstream pope ERROR!\n"; return 1;}
      }
    else if(in_or_out == "try")
      {
      if(verbose) std::cout << "Opening fstream " + p + "\n";
      file.open(path.c_str(), std::ios::in);
      if( !file.is_open() )
        {
        file.close();
        // Checking path to the file exists
        for(pos=0; pos != -1 ; ) {prepos = pos; pos++; pos = path.find('/', pos);}
        sys = "if test -d \"" + path.substr(0, prepos) + "\"; then echo \"T\"; else echo \"F\"; fi";
        pope = popen(sys.c_str(), "r");
        if(fgets(c, 255, pope) == NULL) { std::cout << "Fgets ERROR! Terminating.\n"; std::cin.get(); return 1; }

        if( c[0] == 'T')
          {
          std::cout << "Path not found:" + path + ", creating it.\n";
          file.open(path.c_str(), std::ios::out);
          return 1;
          }
        else if( c[0] == 'F')
          {
          std::cout << "Path not found: " + path + ". Terminating!";
          while( !std::cin.get() ) ;  // ja que sino es pensa que es l'altre cas...
          return 1;
          }
        else {std::cout << "class_fstream pope ERROR!\n"; return 1;}
        }
      }
    else
      { std::cout << "class_fstream default ERROR!\n"; while( !std::cin.get() ) ; return 1; }

    return 0;
    }

  void close() { file.close(); }
};

