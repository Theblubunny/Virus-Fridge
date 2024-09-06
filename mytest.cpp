// Brute force testing script

#include "vdetect.h"
#include <random>
#include <vector>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};


unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);


class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};

class Tester{
public:

    bool insertTestNC(){//testing insert for non-colliding data points

        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        //inserting 4 viruses with known index and non-colliding key & id.
        Virus dataObj = Virus("AATTA", 4987);//index: 5
        result = result && vdetect.insert(dataObj);//testing if insert is successful
        dataObj = Virus("TATAT",5486);//index: 8
        result = result && vdetect.insert(dataObj);
        dataObj = Virus("TTCAT",7739);//index: 21
        result = result && vdetect.insert(dataObj);
        dataObj = Virus("TTGAT",6702);//index: 34
        result = result && vdetect.insert(dataObj);

        result = result && (vdetect.m_currentTable[5].getID() == 4987);
        result = result && (vdetect.m_currentTable[8].getID() == 5486);
        result = result && (vdetect.m_currentTable[21].getID() == 7739);
        result = result && (vdetect.m_currentTable[34].getID() == 6702);

        result = result && vdetect.m_currentSize == 4;

        return result;
    }

    bool insertTestC(){//testing insert for colliding data points & duplicate insert

        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        //inserting 2 viruses with colliding/exact key values. (different IDs)
        //inserting 2 viruses with colliding/exact ID values.  (different Keys)
        Virus dataObj = Virus("AATTA", 4987);//index: 5
        result = result && vdetect.insert(dataObj);//testing if insert is successful
        dataObj = Virus("AATTA",5486);//index: 11
        result = result && vdetect.insert(dataObj);
        dataObj = Virus("TTCAT",7739);//index: 21
        result = result && vdetect.insert(dataObj);
        dataObj = Virus("TTGAT",7739);//index: 34
        result = result && vdetect.insert(dataObj);
        result = result && !vdetect.insert(dataObj);//attempting to insert duplicate virus, expecting to fail.

        //testing if viruses were inserted in the expected buckets
        result = result && (vdetect.m_currentTable[5].getID() == 4987);
        result = result && (vdetect.m_currentTable[11].getID() == 5486);
        result = result && (vdetect.m_currentTable[21].getID() == 7739);
        result = result && (vdetect.m_currentTable[34].getID() == 7739);
        //testing if size changed correctly
        result = result && vdetect.m_currentSize == 4;

        return result;
    }

    bool getVirusTest(){

        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        //inserting non-colliding viruses
        Virus dataObj_1 = Virus("AATTA", 4987);
        result = result && vdetect.insert(dataObj_1);
        Virus dataObj_2 = Virus("TATAT",5486);
        result = result && vdetect.insert(dataObj_2);
        Virus dataObj_3 = Virus("TTCAT",7739);
        result = result && vdetect.insert(dataObj_3);
        Virus dataObj_4 = Virus("TTGAT",6702);
        result = result && vdetect.insert(dataObj_4);

        //test if getVirus() returns the correct virus.
        result = result && (vdetect.getVirus("AATTA", 4987) == dataObj_1);
        result = result && (vdetect.getVirus("TATAT", 5486) == dataObj_2);
        result = result && (vdetect.getVirus("TTCAT", 7739) == dataObj_3);
        result = result && (vdetect.getVirus("TTGAT", 6702) == dataObj_4);

        return result;
    }


    bool getVirusError(){
        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        //inserting 2 non-colliding viruses
        Virus dataObj = Virus("AATTA", 4987);//index: 5
        vdetect.insert(dataObj);//testing if insert is successful
        dataObj = Virus("TATAT",5486);//index: 8
        vdetect.insert(dataObj);

        //test if search returns empty virus if virus does not exist.
        result = result && (vdetect.getVirus("AGATA", 5702) == EMPTY);

        return result;
    }

    bool getVirusColllisionTest(){

        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
        vector<Virus> dataList;
        Random RndID(MINID,MAXID);

        //inserting colliding data
        for (int i=0;i<52;i++){//enough data points to trigger a rehash
            // generating colliding data (same key)
            Virus dataObj = Virus("AATTA", RndID.getRandNum());
            // saving data.
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            vdetect.insert(dataObj);
        }

        //testing if getVirus triggers rehasing.
        for (vector<Virus>::iterator it = dataList.begin(); it != dataList.end(); it++){
            result = result && (*it == vdetect.getVirus((*it).getKey(), (*it).getID()));
            result = result && (vdetect.m_oldTable != nullptr);//if getVirus triggers rehashing
                                                                    // the old table would be deleted.
        }

        return result;
    }

    bool removeTestNC(){
        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        //inserting 4 viruses with known index and non-colliding key & id.
        Virus dataObj_1 = Virus("AATTA", 4987);//index: 5
        vdetect.insert(dataObj_1);//testing if insert is successful
        Virus dataObj_2 = Virus("TATAT",5486);//index: 8
        vdetect.insert(dataObj_2);
        Virus dataObj_3 = Virus("TTCAT",7739);//index: 21
        vdetect.insert(dataObj_3);
        Virus dataObj_4 = Virus("TTGAT",6702);//index: 34
        vdetect.insert(dataObj_4);

        //testing remove
        result = result && vdetect.remove(dataObj_2);
        result = result && vdetect.remove(dataObj_3);
        //testing if size changed correctly
        result = result && (vdetect.m_currNumDeleted == 2);

        //testing if viruses are not found.
        result = result && (vdetect.getVirus(dataObj_2.getKey(), dataObj_2.getID()) == EMPTY);
        result = result && (vdetect.getVirus(dataObj_3.getKey(), dataObj_3.getID()) == EMPTY);

        return result;
    }

    bool removeTestC(){
        bool result = true;
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        //inserting 4 colliding viruses
        Virus dataObj_1 = Virus("AATTA", 4987);//index: 5
        vdetect.insert(dataObj_1);//testing if insert is successful
        Virus dataObj_2 = Virus("AATTA",5486);//index: 11
        vdetect.insert(dataObj_2);
        Virus dataObj_3 = Virus("TTCAT",7739);//index: 21
        vdetect.insert(dataObj_3);
        Virus dataObj_4 = Virus("TTGAT",7739);//index: 34
        vdetect.insert(dataObj_4);


        //testing remove
        result = result && vdetect.remove(dataObj_2);
        result = result && vdetect.remove(dataObj_4);
        //testing if size changed correctly
        result = result && (vdetect.m_currNumDeleted == 2);

        //testing if viruses are not found.
        result = result && (vdetect.getVirus(dataObj_2.getKey(), dataObj_2.getID()) == EMPTY);
        result = result && (vdetect.getVirus(dataObj_4.getKey(), dataObj_4.getID()) == EMPTY);
        //testing if rehash was triggered.
        result = result && vdetect.m_oldTable == nullptr;

        return result;
    }

    bool rehashInsertTriggerTest(){

        Random RndID(MINID,MAXID);
        VDetect vdetect(201, hashCode, DOUBLEHASH);
        bool result = true;
        int insertSize = (vdetect.m_currentCap / 2) + 1;//number of viruses we need to insert to trigger rehash

        //inserting viruses
        for (int i=0;i<insertSize;i++){
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
            // inserting data in to the VDetect object
            vdetect.insert(dataObj);
        }

        //checking if rehashing was triggered
        result = result && vdetect.m_oldTable != nullptr;

        return result;

    }

    bool rehashInsertCompletionTest(){
        Random RndID(MINID,MAXID);
        int tableSize = MINPRIME;
        VDetect vdetect(tableSize, hashCode, DOUBLEHASH);
        bool result = true;

        //number of viruses we need to insert to trigger and complete rehash
        int insertSize = (vdetect.m_currentCap / 2) + 5;

        //inserting viruses
        for (int i=0;i<insertSize;i++){
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum());
            // inserting data in to the VDetect object
            vdetect.insert(dataObj);
        }

        //checking if rehashing was completed
        result = result && vdetect.m_oldTable == nullptr;
        //checking if curr table has the current size
        result = result && vdetect.m_currentCap == vdetect.findNextPrime(tableSize * 4) ;;
        //checking if all viruses were transferred properly
        result = result && vdetect.m_currentSize == insertSize;


        return result;
    }

    bool rehashRemoveTriggerTest(){

        Random RndID(MINID,MAXID);
        int tableSize = 101;
        VDetect vdetect(tableSize, hashCode, DOUBLEHASH);
        bool result = true;

        //inserting viruses
        Virus dataObj_1 = Virus("AATTA", 4987);//index: 5
        vdetect.insert(dataObj_1);//testing if insert is successful
        Virus dataObj_2 = Virus("TATAT",5486);//index: 8
        vdetect.insert(dataObj_2);
        Virus dataObj_3 = Virus("TTCAT",7739);//index: 21
        vdetect.insert(dataObj_3);
        Virus dataObj_4 = Virus("TTGAT",6702);//index: 34
        vdetect.insert(dataObj_4);
        Virus dataObj_5 = Virus("CAAAA",9281);//index: 57
        vdetect.insert(dataObj_5);
        Virus dataObj_6 = Virus("CTACG",3654);//index: 71
        vdetect.insert(dataObj_6);

        //removing viruses
        result = result && vdetect.remove(dataObj_1);
        result = result && vdetect.remove(dataObj_2);
        result = result && vdetect.remove(dataObj_3);
        result = result && vdetect.remove(dataObj_4);
        result = result && vdetect.remove(dataObj_5);

        //checking if rehashing was triggered
        result = result && vdetect.m_currentCap == vdetect.findNextPrime(tableSize * 4);

        return result;

    }

    bool rehashRemoveCompletionTest(){

        Random RndID(MINID,MAXID);
        int tableSize = 101;
        VDetect vdetect(tableSize, hashCode, DOUBLEHASH);
        bool result = true;

        //inserting viruses
        Virus dataObj_1 = Virus("AATTA", 4987);//index: 5
        vdetect.insert(dataObj_1);//testing if insert is successful
        Virus dataObj_2 = Virus("TATAT",5486);//index: 8
        vdetect.insert(dataObj_2);
        Virus dataObj_3 = Virus("TTCAT",7739);//index: 21
        vdetect.insert(dataObj_3);
        Virus dataObj_4 = Virus("TTGAT",6702);//index: 34
        vdetect.insert(dataObj_4);
        Virus dataObj_5 = Virus("CAAAA",9281);//index: 57
        vdetect.insert(dataObj_5);
        Virus dataObj_6 = Virus("CTACG",3654);//index: 71
        vdetect.insert(dataObj_6);

        //removing viruses
        result = result && vdetect.remove(dataObj_1);
        result = result && vdetect.remove(dataObj_2);
        result = result && vdetect.remove(dataObj_3);
        result = result && vdetect.remove(dataObj_4);
        result = result && vdetect.remove(dataObj_5);

        //checking if rehashing was triggered
        result = result && vdetect.m_currentCap == vdetect.findNextPrime(tableSize * 4);
        //checking if old table has been removed
        result = result && vdetect.m_oldTable == nullptr;

        return result;

    }



};

int main(){

    Tester tester;

    if(tester.insertTestNC()){
        cout << "PASSED - Non-colliding data points insert test." << endl;
    }else{
        cout << "FAILED - Non-colliding data points insert test." << endl;
    }
    if(tester.insertTestC()){
        cout << "PASSED - Colliding data points & duplicate insert test." << endl;
    }else{
        cout << "FAILED - Colliding data points & duplicate insert test." << endl;
    }

    if(tester.getVirusTest()){
        cout << "PASSED - getVirus() on few non-colliding keys." << endl;
    }else{
        cout << "FAILED - getVirus() on few non-colliding keys." << endl;
    }
    if(tester.getVirusError()){
        cout << "PASSED - getVirus() on non-existent virus." << endl;
    }else{
        cout << "FAILED - getVirus() on non-existent virus." << endl;
    }
    if(tester.getVirusColllisionTest()){
        cout << "PASSED - getVirus() does not trigger rehash on colliding data points." << endl;
    }else{
        cout << "FAILED - getVirus() triggers rehash on colliding data points." << endl;
    }

    if(tester.removeTestNC()){
        cout << "PASSED - Non-colliding data points remove test." << endl;
    }else{
        cout << "FAILED - Non-colliding data points remove test." << endl;
    }
    if(tester.removeTestC()){
        cout << "PASSED - Colliding data points remove test." << endl;
    }else{
        cout << "FAILED - Colliding data points remove test." << endl;
    }

    if(tester.rehashInsertTriggerTest()){
        cout << "PASSED - Insert rehash was triggered." << endl;
    }else{
        cout << "FAILED - Insert rehash was not triggered." << endl;
    }
    if(tester.rehashInsertCompletionTest()){
        cout << "PASSED - Insert rehash was completed." << endl;
    }else{
        cout << "FAILED - Insert rehash was not completed." << endl;
    }
    if(tester.rehashRemoveTriggerTest()){
        cout << "PASSED - Remove rehash was triggered." << endl;
    }else{
        cout << "FAILED - Remove rehash was not triggered." << endl;
    }
    if(tester.rehashRemoveCompletionTest()){
        cout << "PASSED - Remove rehash was completed." << endl;
    }else{
        cout << "FAILED - Remove rehash was not completed." << endl;
    }

    return 0;
}

unsigned int hashCode(const string str) {
    unsigned int val = 0 ;
    const unsigned int thirtyThree = 33 ;  // magic number from textbook
    for ( int i = 0 ; i < str.length(); i++)
        val = val * thirtyThree + str[i] ;
    return val ;
}

string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}
