/*****************************************************************************
 *   GATB software relying on Genome Analysis Toolbox with de-Bruijn graph 
 *   Copyright (C) 2014-2016  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _TOOL_ItrAs_HPP_
#define _TOOL_ItrAs_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <gatb/debruijn/impl/GraphUnitigs.hpp>
#include <gatb/debruijn/impl/Simplifications.hpp>
#include <cstring>
#include <vector>
#include <thread>
#include <sys/sysinfo.h>
#include <dirent.h>
#include <iomanip>
#include <iostream>
#include <fstream>
/********************************************************************************/
static const uint32_t kMaxLine = (1 << 20);
static char buffer[kMaxLine];
static std::string FormatString(const char *fmt, ...){
  va_list ap;
  va_start(ap, fmt);
  vsprintf(buffer, fmt, ap);
  va_end(ap);
  return buffer;
}

static void MakeDir(const std::string &directory){
  DIR *dir = opendir(directory.c_str());
  if (dir == NULL){
    if (mkdir(directory.c_str(), 0777) == -1){
      std::cout<<"can't create directory\n";
      exit(1);
    }
  }
}

static void RemoveTmpFiles(const std::string &directory, const std::string &oldContigFile){
  ///Remove *.fa.glue* files
  std::string cmd = "rm -f "+directory+"/*.fa.glue* "+
                    directory+"/*.unitigs.fa "+directory+"/*.h5 " + oldContigFile;
  int err = std::system(cmd.c_str());
  if(err) ;//So the compiler doesn't complain at us.
}

static void FinalizeContigFile(const std::string &directory, const std::string &contigFile){
  std::string cmd = "mv "+contigFile+ " "+directory+"/contig.fa";
  int err = std::system(cmd.c_str());
  if(err) ;//So the compiler doesn't complain at us.
}
/********************************************************************************/

class ItrAs : public Tool
{
	bool keepIsolatedTigs       ;
  u_int64_t nbSmallContigs    ;
  u_int64_t totalNt           ;
  u_int64_t maxContigLen      ;
  template <typename Graph_type, typename Node, typename Edge, size_t span>
  void assembleFrom(Node startingNode, Graph_type& graph, IBank *outputBank);
public:
/****************Vars Used for Options******************/
	///used for karect
  std::string celltype;
  std::string matchtype;
  std::string aggressive;
  std::string precor_stages;
  std::string precor_trim;
  bool bool_precorrect;
 	///used for assembly
  int mink;
  int maxk;
  int step;
  unsigned int min_contig;
  std::string directory;
  std::string read_file;
/**************End of vars used for Options*************/
  u_int64_t nbContigs;
  std::vector<int> kVals;

	// Constructor
	ItrAs ();

	// Actual job done by the tool is here
	void execute ();

	//puts all the options into our variables
	void ParseOptions();

	//set the kvals used for assembly
	void SetKVals(){ kVals.clear(); for(int i = mink; i <= maxk; i += step){ kVals.push_back(i); } }

	//Used to do the actual assembly (Minia)
	template <typename Graph_type, typename Node, typename Edge, size_t span>
  string assemble (Graph_type& graph);
};

/********************************************************************************/

#endif /* _TOOL_ItrAs_HPP_ */

