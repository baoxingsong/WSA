//
// Created by Baoxing Song on 2019-04-02.
//

#include "fasta.h"

void charByCharRead(std::istream& stream, std::map<std::string, std::string> & seqs ) {  // not using, could be deleted
//    if (!stream.good()) throw EggRuntimeError("cannot import fasta data: invalid stream");

    // ignores everything until the first >
    char c=' ';
    while (c!='>') {
        c= stream.get();
        if (stream.eof()) {
            // accepting empty files
            return;
        }
    }

    // main loop over sequences (each starts after having read a >)
    while (stream.good()) {

        std::string name;
        std::string sequence;

        // gets the name
        c= stream.get();
        while (c!='\n') {
            if (c!='\r') name.push_back(c);
            c= stream.get();
        }

        name.push_back('\0');

        // gets the sequence
        c = stream.get();
        while (c!='>' && !stream.eof()) {
            char d = (c!='\n' && c!='\r' && c!='\r')?c:'\0';
            if (d) sequence.push_back(d);
            c= stream.get();
        }
        sequence.push_back('\0');
        seqs[name]=sequence;
    }
}


void readFastaFile( const std::string& filePath, std::map<std::string, std::string>& sequences ){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }
    std::string name="";
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){ // changed on 19 Sep, 2018, by reading the source code from Heng, Li.
            // this method is faster than regular expression
            if( name.size()>0 ){
                std::string sequence = sequencestream.str();
                sequences[name]=sequence;
            }
            name=line.substr(1, line.find(" ", 0)-1);
            if( name[0] == '>' ){
                name=name.substr(1, name.find("\t", 0)-1);
            }else{
                name=name.substr(0, name.find("\t", 0)-1);
            }
            sequencestream.str(std::string());
        }else{
            sequencestream << line;
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        if( sequence.size()>0 ){
            sequences[name]=sequence;
        }
    }
    infile.close();
}


void charByCharCategoryRead( const std::string& filePath, std::map<std::string, std::string> & seqs, std::map<std::string, int16_t * > & weight ) {
    for( std::map<std::string, std::string>::iterator it = seqs.begin(); it!=seqs.end(); ++it ){
        weight[it->first] = new int16_t[it->second.size()];
    }
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }
    std::string name="";
    std::string line;
    int16_t thisCategory;
    int32_t thisId;
    while (std::getline(infile, line)){
        if( line[0] == '>'  ){
            name=line.substr(1, line.find(" ", 0)-1);
            if( name[0] == '>' ){
                name=name.substr(1, name.find("\t", 0)-1);
            }else{
                name=name.substr(0, name.find("\t", 0)-1);
            }
            thisId=0;
        }else{
            std::stringstream ss(line);
            std::string item;
            while (std::getline(ss, item, '\t')) {
                if (item.length() > 0) {
                    weight[name][thisId] = std::stoi(item);
                    ++thisId;
                }
            }
        }
    }
}
