/*
 *  class_fstream.h
 *
 *  Aquest header conte class_fstream, que gestiona els arxius d'entrada i sortida awsomelyment
 *
 */

#include "includes.h"

class class_fstream
{  // Conte el fstream i la funcio per obrirlo be (amb path = "/home/ignasi/IDIBAPS/Alignments/"+p)
public:
  std::string path, in_or_out, tag;
  std::fstream file;

  bool obre(std::string p, std::string in_or_out, std::string t = "", int verbose=1)
    {  // Pero s'obre nomes donantli el path corresponent
    int pos, prepos;
    char c[255];
    std::string sys;
    FILE *pope;

    path = p;

    if( t == "") tag = p;
    else tag = t;

    if(in_or_out == "in")
      {
      if(verbose) std::cout << "Obrint input " + tag + "\n";
      file.open(path.c_str(), std::ios::in);
      if( !file.is_open() )
        {
        file.close();
        std::cout << "Impossible obrir " + path + ". Terminating!";
        return 1;
        }
      }
    else if(in_or_out == "out")
      {
      if(verbose) std::cout << "Obrint output " + p + "\n";

      // Comprovem que existeixi el path a l'arxiu
      for(pos=0; pos != -1 ; ) {prepos = pos; pos++; pos = path.find('/', pos);}
      sys = "if test -d \"" + path.substr(0, prepos) + "\"; then echo \"T\"; else echo \"F\"; fi";
      pope = popen(sys.c_str(), "r"); if(fgets(c, 255, pope) == NULL) std::cout << "Fgets ERROR!\n";

      if( c[0] == 'T') file.open(path.c_str(), std::ios::out);
      else if( c[0] == 'F')
        {
        std::cout << "No existeix el path a " + path + ". Terminating!";
        return 1;
        }
      else {std::cout << "class_fstream pope fucked up!\n"; return 1;}
      }
    else if(in_or_out == "try")
      {
      if(verbose) std::cout << "Obrint fstream " + p + "\n";
      file.open(path.c_str(), std::ios::in);
      if( !file.is_open() )
        {
        file.close();
        // Comprovem que existeixi el path a l'arxiu
        for(pos=0; pos != -1 ; ) {prepos = pos; pos++; pos = path.find('/', pos);}
        sys = "if test -d \"" + path.substr(0, prepos) + "\"; then echo \"T\"; else echo \"F\"; fi";
        pope = popen(sys.c_str(), "r");
        if(fgets(c, 255, pope) == NULL) { std::cout << "Fgets ERROR! Terminating.\n"; std::cin.get(); return 1; }

        if( c[0] == 'T')
          {
          std::cout << "No hem trobat " + path + ". L'haurem de construir.\n";
          file.open(path.c_str(), std::ios::out);
          return 1;
          }
        else if( c[0] == 'F')
          {
          std::cout << "No existeix el path a " + path + ". Terminating!";
          while( !std::cin.get() ) ;  // ja que sino es pensa que es l'altre cas...
          return 1;
          }
        else {std::cout << "class_fstream pope fucked up!\n"; return 1;}
        }
      }
    else
      { std::cout << "class_fstream diu: burru!\n"; while( !std::cin.get() ) ; return 1; }

    return 0;
    }

  void tanca() { file.close(); }
};

