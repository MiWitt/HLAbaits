/* 
 * File:   main.cpp
 * Author: mwittig
 *
 * Created on May 21, 2013, 8:18 AM
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <string.h>


#include "CMyException.h"
#include "CMultiFasta.h"

using namespace std;

#define INCLUDE_FLANKING_BASES 200
// #define PROBE_SIZE 120
int PROBE_SIZE = 120;
#define MIN_NEW_BAIT_SIZE 30
// #define MISMATCH 4
int MISMATCH = 5;
int TILING_STEP_SIZE = 120;

string SCRATCH_DIR =".";

// alle blöcke die größer MAX_EXON_SIZE sind behalten nur ihre 500bp flanken
// plotten von v3.12 gab maximale exongröße 599, somit sind wir hier sicher im Intron
// #define MAX_EXON_SIZE 650
#define MAX_EXON_SIZE 500650
#define INCLUDE_INTRON_FLANKS 150



/*
 * 
 */


vector<string> Tokenize(const string& str, const string& delimiters = "\t")
{
    vector<string> tokens;
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        if(str.length() > pos && delimiters.find(str[pos]) != string::npos )
            lastPos=pos+1;
        else
            lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
}


string cmd_exec(const string& cmd)
{
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        throw("ERROR from popen()");

    char buffer[128];
    ostringstream result("");

    while(!feof(pipe))
    {
        if(fgets(buffer, 128, pipe) != NULL)
           result << buffer;
    }
    pclose(pipe);

    return result.str();
}

vector<pair<int,int> > parseResults(const string& filename)
{
    try
    {
        vector<pair<int,int> > vctReturn;
        ifstream in(filename.c_str());
        if(!in)
        {
            cerr << "Can't open " << filename << endl;
             exit(EXIT_FAILURE);
        }
        size_t truncated_reads = 0;
        string strLine;
        while(getline(in,strLine))
            if(strLine.length() > 0 && strLine[0] == '-')
                break;
        while(getline(in,strLine))
        {
            vector<string> parsed = Tokenize(strLine);
            //  0      1       2       3        4       5       6       7       8     9        10     11      12      13      14      15      16      17      18         19              20       21
            //match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q        Q       Q       Q       T       T       T       T     block   blockSizes  qStarts     tStarts
            //match   match           count   bases   count   bases                  name    size    start   end     name    size    start   end     count
            int RefEnd      = atoi(parsed[16].c_str())-1;
            int RefStart    = atoi(parsed[15].c_str());
            vctReturn.push_back(pair<int,int>(RefStart,RefEnd));
        }
        in.close();

        return vctReturn;
    }
    catch(...)
    {
        cerr << "Error parsing psl" << endl;
        exit(EXIT_FAILURE);
    }
}

vector<pair<int,int> >  getUnmappedRegions(const string& matches)
{
    vector<pair<int,int> >  vctRet;
    size_t starter = matches.find_first_of("0",0);

    while(starter != string::npos)
    {
        int i = starter;
        int counter =0;
        for(; i < matches.size() && matches[i] == '0';i++)
            counter++;
        if(counter >= MIN_NEW_BAIT_SIZE)
            vctRet.push_back(pair<int,int>((int)starter,(int)starter+counter-1));
        starter = matches.find_first_of("0",i);
    }
    return vctRet;
}

map<string,string>  tile(string ID, string seq, int start, int end, int step_size = TILING_STEP_SIZE)
{
    map<string,string> mapRet;
    cerr << "tiling: " << setw(20) << right << ID << ", uncovered bases:" << setw(10) << right << end-start+1 << endl;
    while( (end-start)+1 < PROBE_SIZE)
    {
        if(rand() % 2 == 0)
        {
            if(start > 0)
                start--;
            if((end-start)+1 < PROBE_SIZE && end < seq.size()-1)
                end++;
        }
        else
        {
            if(end < seq.size()-1)
                end++;
            if((end-start)+1 < PROBE_SIZE && start > 0)
                start--;

        }
    }
    if((end-start)+1 == PROBE_SIZE)
    {
        ostringstream osr("");
        osr << ID << '_' << start << '_' << end;
        mapRet[osr.str()]=seq.substr(start,(end-start)+1);
    }
    else
    {
        if(rand() % 2 == 1)
        {
            ostringstream osr("");
            osr << ID << '_' << start << '_' << start+PROBE_SIZE-1;
            mapRet[osr.str()]=seq.substr(start,PROBE_SIZE);
            map<string,string> mapTmp = tile(ID,seq,start+step_size,end);
            mapRet.insert(mapTmp.begin(),mapTmp.end());
        }
        else
        {
            ostringstream osr("");
            osr << ID << '_' << (end-PROBE_SIZE)+1 << '_' << end;
            mapRet[osr.str()]=seq.substr((end-PROBE_SIZE)+1,PROBE_SIZE);
            map<string,string> mapTmp = tile(ID,seq,start,end-step_size);
            mapRet.insert(mapTmp.begin(),mapTmp.end());
        }
    }
    return mapRet;
}

int main(int argc, char** argv) {

    try
    {
        // Achtung Anpassung auch weiter unten: if(argc==7) und ifstream inProbes(argv[6]);
        if(argc != 6 && argc != 7)
        {
            cout << "run like: " << basename(argv[0]) << " probe_length mis_matches tiling_step ref_seq_fasta scratch [existing_probes_fasta]" << endl;
            cout << "probe_length: define the length of the probes here probes (recommended 120)"  
                 << "mis_matches: for each new fasta sequence all existing probes are aligned to it using blat. Define the number of allowed mismatches here (recommended 5)"  
                 << "tiling_step: for nun overlapping probes set tiling step = probe_length. Set it to probe_length/2 for half overlapping (recommended 60)"  
                 << "ref_seq_fasta: Provide your collection of fasta sequences here. Set the path to a multiple fasta file "  
                 << "scratch: Set the path to a scratch directory where " << basename(argv[0]) << " has read write access "   
                 << "[existing_probes_fasta]: Provide a file of already existing probes here. Be aware of the important note below"  
                 << "All output goes to stdout, so redirect the output to your destination file. (e.g. > /foo/bar/my_baits.fa)" << endl;
            cout << "IMPORTANT: if you provide an existing probes (baits) file, provide it in following format. First line unique probe name, second line probe sequence, third line unique probe name, etc. NO BLANK LINES, NO REDUNDANT PROBE NAMES." << endl;
            exit(EXIT_FAILURE);
        }
        PROBE_SIZE = atoi(argv[1]);
        MISMATCH = atoi(argv[2]);
        TILING_STEP_SIZE = atoi(argv[3]);
        SCRATCH_DIR = string(argv[5])+"/";

        string scratch_template = SCRATCH_DIR+"template_seq.fasta";
        string scratch_baits = SCRATCH_DIR+"template_baits.fasta";
        string scratch_psl = SCRATCH_DIR+"output.psl";
        string ref_seq_fasta = argv[4];

        string strLine;
        string strLine2;

        map<string,string> new_baits;
        map<string,string> newer_baits;
        string ID = "5884b7a6-2879-4e85-841b-9284beb34974";
        vector<string>  sequences;

        if(argc==7)
        {
            ifstream inProbes(argv[6]);
            while(getline(inProbes,strLine) && getline(inProbes,strLine2))
            {
                new_baits[strLine]=strLine2;
            }
            inProbes.close();
        }
        
        
        CMultiFasta mfa(ref_seq_fasta.c_str());
        for(CMultiFasta::iterator iMfa = mfa.begin(); iMfa != mfa.end(); iMfa++)
        {
            strLine = iMfa->first;
            strLine2 = iMfa->second;
            
            vector<string> tok = Tokenize(strLine," ");
            if(tok[0] != ID)
            {
                string tiling_seq="";
                for(int i = 0; i < sequences.size(); i++)
                {
                    if(sequences[i].size() > MAX_EXON_SIZE)
                        sequences[i].erase(INCLUDE_INTRON_FLANKS,sequences[i].size()-(2*INCLUDE_INTRON_FLANKS));
                    tiling_seq+=sequences[i];
                }
                ofstream out(scratch_template.c_str());
                out << ID << endl;
                out << tiling_seq << endl;
                out.close();

                out.open(scratch_baits.c_str());
                for(map<string,string>::iterator iter = new_baits.begin();iter != new_baits.end();iter++)
                {
                    out << iter->first << endl << iter->second << endl;
                }
                out.close();

                sequences.clear();
                sequences.push_back(strLine2);

                ostringstream osr("");
                osr << "blat  -minScore=" << PROBE_SIZE-MISMATCH 
                    <<  " -repMatch=64000 -fine -fastMap " << scratch_template
                    << " " << scratch_baits << " " << scratch_psl;
                //cout << osr.str() << endl;
                cmd_exec(osr.str());
                // blat -minScore=77 -repMatch=64000 -fine -fastMap HLA_v3.12_block_wise.fasta std_referenz_baits.txt stdout

                vector<pair<int,int> > mapping_res = parseResults(scratch_psl);
                string matches(tiling_seq.size(),'0');
                for(int i = 0; i < mapping_res.size();i++)
                    for(int j = mapping_res[i].first; j < mapping_res[i].second;j++)
                        matches[j]='1';

                vector<pair<int,int> > unmapped = getUnmappedRegions(matches);

                for(int i = 0; i < unmapped.size();i++ )
                {
                    if(tiling_seq.size() < PROBE_SIZE)
                        continue;
                    map<string,string> add_baits = tile(ID,tiling_seq,unmapped[i].first,unmapped[i].second,TILING_STEP_SIZE);
                    new_baits.insert(add_baits.begin(),add_baits.end());
                    if(argc==6)
                        newer_baits.insert(add_baits.begin(),add_baits.end());
                }

                ID=tok[0];

                //for(int i = 0; i < unmapped.size();i++ )
                //    cout << ID << " " << unmapped[i].first << " - " << unmapped[i].second << endl;
                cerr << new_baits.size() << " baits after " << ID << endl;
            }
            else
                sequences.push_back(strLine2);
        }
        


        if(argc==6)
            for(map<string,string>::iterator iter = newer_baits.begin(); iter != newer_baits.end(); iter++)
                cout << iter->first << endl << iter->second << endl;
        else
            for(map<string,string>::iterator iter = new_baits.begin(); iter != new_baits.end(); iter++)
                cout << iter->first << endl << iter->second << endl;
    
        return (EXIT_SUCCESS);
    }
    catch(CMyException e)
    {
        cerr << "An error occured: " << e.what() << endl;
    }
    catch(exception e)
    {
        cerr << "An error occured: " << e.what() << endl;
    }
    catch(...)
    {
        cerr << "An undefined error occured ..." << endl;
    }
    return  EXIT_FAILURE;
}

