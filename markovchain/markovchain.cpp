 /* Order-2 Markov chain program
    by Andrew M. Launder
    
    main() loop adapted from a Python example
    (see http://code.activestate.com/recipes/194364-the-markov-chain-algorithm/ ).
    
    See README for proper code usage. */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::string> LineSplit(std::ifstream& text) {
 // Splits a string and stores it in a vector.
    std::string line, word;
    std::getline(text, line);
    std::istringstream linestream(line);
    std::vector<std::string> linesplit;
    
    while (linestream >> word) {
        linesplit.push_back(word);
    }
    
    return linesplit;
}

unsigned long int NLines(char* filename) {
 // Determine number of lines in input text.
    std::ifstream text(filename);
    
    if (!text) {
        std::cerr << "No input file found under that name!" << std::endl;
        return 0;
    }
    
    unsigned long int nlines = 0;
    
    std::string line;
    
    while (std::getline(text, line)) {
        ++nlines;
    }
    
    return nlines;
}

unsigned long int NWords(char* filename, unsigned long int nlines) {
 // Determine number of words in input text.
    std::ifstream text(filename);
    
    if (!text) {
        return 0;
    }
     
    unsigned long int words, nwords;
    nwords = 0;
    
    std::vector<std::string>* linesplit = new std::vector<std::string>;
    
    unsigned long int i;
    for (i = 0; i < nlines; ++i) {
        *linesplit = LineSplit(text);
        words = linesplit->size();
        nwords += words;
    }
    delete linesplit;
    
    return nwords;
}

std::vector<std::vector<std::string>> MarkovData(char* filename, unsigned long int nlines, unsigned long int nwords) {
 /* Parse data into vector of vectors containing every iteration of
    three consecutive words. Make sure all empty lines are removed from
    input text. */
    std::ifstream text(filename);
    
    if (!text) {
        std::vector<std::string> errorvector(1);
        errorvector[0] = "Error string.";
        std::vector<decltype(errorvector)> errorreturn(1);
        errorreturn[0] = errorvector;
        return errorreturn;
    }
    
    std::string w1, w2, w3;
    w1 = w2 = w3 = "\n";
    unsigned int words, totalwords;
    totalwords = 0;
    std::vector<std::string>* ws = new std::vector<std::string>(3);
    std::vector<std::string>* linesplit = new std::vector<std::string>;
    std::vector<std::vector<std::string>> markovdata(nwords);
    
    unsigned long int i, j;
    for (i = 0; i < nlines; ++i) {
        *linesplit = LineSplit(text);
        words = linesplit->size();
        for (j = 0; j < words; ++j) {
            w3 = (*linesplit)[j];
            *ws = {w1, w2, w3};
            markovdata[totalwords + j] = *ws;
            w1 = w2;
            w2 = w3;
        }
        totalwords += words;
    }
    delete ws;
    delete linesplit;
    
    return markovdata;
}

int main(int argc, char** argv) {
 // Markov chain algorithm.
    unsigned long int maxwords;
    char* filename;
    
    if (argc <= 2) {
        maxwords = 100;
        filename = argv[1];
    }
    else {
        std::istringstream maxwordsstr(argv[1]);
        maxwordsstr >> maxwords;
        if (maxwords <= 0 || maxwords >= 10000) {
            std::cerr << "Word count must be greater than 0 and less than 10000." << std::endl;
            return EXIT_FAILURE;
        }
        
        filename = argv[2];
    }
    
    if (!filename) {
        std::cerr << "No input file given!" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string w1, w2, w3;
    w1 = w2 = "\n";
    std::vector<std::string>* ws = new std::vector<std::string>(3);
    unsigned long int nlines, nwords;
    nlines = NLines(filename);
    nwords = NWords(filename, nlines);
    std::vector<std::vector<std::string>>* markovdata = new std::vector<std::vector<std::string>>;
    *markovdata = MarkovData(filename, nlines, nwords);
    
    if ((*markovdata)[0][0] == "Error string.") {
        return EXIT_FAILURE;
    }
    
    std::cout << std::endl;
    std::cout << maxwords << "-word text generated by order-2 Markov chain: " << std::endl;
    std::cout << std::endl;
    
    std::vector<std::vector<std::string>>* matches = new std::vector<std::vector<std::string>>;
    std::vector<std::string>* match = new std::vector<std::string>(3);
    unsigned long int matchnumber, w3number;
    
    unsigned long int i, j;
    for (i = 0; i < maxwords; ++i) {
        for (j = 0; j < nwords; ++j) {
            *match = (*markovdata)[j];
            if (w1 == (*match)[0] && w2 == (*match)[1]) {
                matches->push_back(*match);
            }
        }
        matchnumber = matches->size();
        srand(time(NULL));
        w3number = rand() % matchnumber;
        *ws = (*matches)[w3number];
        w3 = (*ws)[2];
        std::cout << w3 << " ";
        w1 = w2;
        w2 = w3;
        matches->clear();
    }
    delete ws;
    delete markovdata;
    delete matches;
    delete match;
    std::cout << std::endl;
    std::cout << std::endl;
    
    std::cout << "Program exiting successfully." << std::endl;
    std::cout << std::endl;
    
    return 0;
}
