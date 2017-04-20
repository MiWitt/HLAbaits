/* 
 * File:   CMultiFasta.cpp
 * Author: mwittig
 * 
 * Created on July 20, 2012, 1:00 PM
 */
#include <cstdio>
#include <cstdlib>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "CMyException.h"
#include "CMultiFasta.h"

#define BOOST_ALL_NO_LIB

CMultiFasta::CMultiFasta()
{
   
}

CMultiFasta::CMultiFasta(const string& filename)
{
    readSequence(filename);
   
}

CMultiFasta::CMultiFasta(const CMultiFasta& orig)
{

}

CMultiFasta::~CMultiFasta()
{

}

void CMultiFasta::readSequence(const string& filename)
{
    try
    {
        ifstream in(filename.c_str());
        if(!in)
            throw string("Can't open ")+filename;
        string strLine;
        string strHeader = "";
        ostringstream strSeq("");

        while(getline(in,strLine))
        {
            if(strLine.size() == 0)
                continue;
            if(strLine[0] == '>')
            {

                if(strHeader!="")
                {
                    m_mapSequences[strHeader]=strSeq.str();
                }
                else
                    CMyException("empty header. discard entries and continue with ")+strLine;
                strHeader=strLine.substr(1,strLine.size());
                strSeq.str("");
            }
            else
                strSeq << strLine;
        }
        in.close();
        if(strHeader!="")
        {
            m_mapSequences[strHeader]=strSeq.str();
        }

    }
    catch(const string& err)
    {
        throw CMyException(err);
    }
    catch(const CMyException& err)
    {
        throw err;
    }
    catch(...)
    {
        throw CMyException("unknown error in CMultiFasta::readSequence ");
    }
}


string CMultiFasta::operator[](const string& val)
{
    map<string,string>::iterator iter = m_mapSequences.find(val);
    if(m_mapSequences.end() == iter)
        return "";
    return iter->second;
}


