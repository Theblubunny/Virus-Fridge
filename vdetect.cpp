#include "vdetect.h"
VDetect::VDetect(int size, hash_fn hash, prob_t probing = DEFPOLCY){

    if(size > MAXPRIME) m_currentCap = MAXPRIME;//Initializing size to be prime
    else if(size < MINPRIME) m_currentCap = MINPRIME;
    else if (!isPrime(size)) m_currentCap = findNextPrime(size);
    else m_currentCap = size;

    m_currProbing = probing;//Initializing curr table variables
    m_hash = hash;
    m_currentTable = new Virus[m_currentCap];
    m_currentSize = 0;
    m_currNumDeleted = 0;

    m_oldTable = nullptr;//Initializing old table variables
    m_oldProbing = NONE;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;

    m_newPolicy = NONE;

}

VDetect::~VDetect(){
    delete [] m_currentTable;;
    delete [] m_oldTable;
}

void VDetect::changeProbPolicy(prob_t policy){
    m_newPolicy = policy;
}

bool VDetect::insert(Virus virus){

    if(virus.getID() > MAXID || virus.getID() < MINID  || !(getVirus(virus.getKey(), virus.getID()) == EMPTY)){
        return false;   //Protects against invalid IDs and duplicate inserts.
    }

    int hashValue = m_hash(virus.getKey()) % m_currentCap;
    int counter = 0;
    bool insertBool = false;//bool used to keep track of whether insertion has been completed.

    //if space is not empty
    //if IDs don't match
    //if space is not deleted
    while(!m_currentTable[hashValue].getKey().empty() && m_currentTable[hashValue].getID() != virus.getID()
            && !(m_currentTable[hashValue] == DELETED)){
        if (m_currProbing == QUADRATIC) {
            hashValue = ((hashValue % m_currentCap) + (counter * counter)) % m_currentCap;
            counter++;
        } else if (m_currProbing == DOUBLEHASH) {
            hashValue = ((hashValue % m_currentCap) + counter * (11 - (hashValue % 11))) % m_currentCap;
            counter++;
        }
    }

    //insert if space is empty and not filled by deleted virus.
    if(m_currentTable[hashValue].getKey().empty() || m_currentTable[hashValue] == DELETED){
        m_currentTable[hashValue] = virus;
        m_currentSize++;
        insertBool = true;
    }

    //rehash if rehashing has already begun OR load factor is greater than 50%
    if((insertBool && m_oldTable != nullptr) || (insertBool && lambda() > 0.5)){
        rehash();
    }

    return insertBool;
}

bool VDetect::remove(Virus virus){

    bool removeBool = false;

    string key = virus.getKey();
    int hashValue = m_hash(virus.getKey()) % m_currentCap;
    int counter = 0;

    //checking for virus in curr table
    while (!m_currentTable[hashValue].getKey().empty()
           && !(m_currentTable[hashValue] == virus)) {
        if (m_currProbing == QUADRATIC) {
            hashValue = ((hashValue % m_currentCap) + (counter * counter)) % m_currentCap;
            counter++;
        } else if (m_currProbing == DOUBLEHASH) {
            hashValue = ((hashValue % m_currentCap) + counter * (11 - (hashValue % 11))) % m_currentCap;
            counter++;
        }
    }

    //if found, tag it as DELETED
    if (m_currentTable[hashValue] == virus) {//check if correct key was found
        m_currentTable[hashValue] = DELETED;
        m_currNumDeleted++;
        removeBool = true;
    }

    //if the node wasn't found in current table, check old table if it exists
    if (m_oldTable != nullptr) {
        int oldHash = m_hash(key) % m_oldCap;
        counter = 0;

        //checking old table.
        while (!m_oldTable[oldHash].getKey().empty()
               && !(m_oldTable[oldHash] == virus)) {
            if (m_oldProbing == QUADRATIC) {
                oldHash = ((oldHash % m_oldCap) + (counter * counter)) % m_oldCap;
                counter++;
            } else if (m_oldProbing == DOUBLEHASH) {
                oldHash = ((oldHash % m_oldCap) + counter * (11 - (oldHash % 11))) % m_oldCap;
                counter++;
            }
        }
        if (m_oldTable[oldHash] == virus ) {//check if correct key was found
            m_oldTable[oldHash] = DELETED;
            m_oldNumDeleted++;
            removeBool = true;
        }
    }

    //rehash if rehashing has already begun || deleted ration 80%
    if((removeBool && m_oldTable != nullptr) || (removeBool && deletedRatio() > 0.8)){
        rehash();
    }

    return removeBool;
    
}

Virus VDetect::getVirus(string key, int id) const {

    int hashValue = m_hash(key) % m_currentCap;
    int counter = 0;

    //checking curr table for virus.
    while (!m_currentTable[hashValue].getKey().empty()
           && m_currentTable[hashValue].getID() != id) {
        if (m_currProbing == QUADRATIC) {
            hashValue = ((hashValue % m_currentCap) + (counter * counter)) % m_currentCap;
            counter++;
        } else if (m_currProbing == DOUBLEHASH) {
            hashValue = ((hashValue % m_currentCap) + counter * (11 - (hashValue % 11))) % m_currentCap;
            counter++;
        }
    }

    //if virus is found, return it.
    if (m_currentTable[hashValue].m_key == key && m_currentTable[hashValue].m_id == id) {//check if correct key was found
        return m_currentTable[hashValue];
    }


    //if the node wasn't found in current table, check old table if it exists
    if (m_oldTable != nullptr) {
        int oldHash = m_hash(key) % m_oldCap;
        counter = 0;

        while (!m_oldTable[oldHash].getKey().empty()
               && m_oldTable[oldHash].getID() != id) {
            if (m_oldProbing == QUADRATIC) {
                oldHash = ((oldHash % m_oldCap) + (counter * counter)) % m_oldCap;
                counter++;
            } else if (m_oldProbing == DOUBLEHASH) {
                oldHash = ((oldHash % m_oldCap) + counter * (11 - (oldHash % 11))) % m_oldCap;
                counter++;
            }
        }

        if (m_oldTable[oldHash].getKey() == key && m_oldTable[oldHash].getID() == id ) {//check if correct key was found
            return m_oldTable[oldHash];
        }
    }

    //if virus is not found, create empty virus and return it
    Virus emptyVirus = EMPTY;
    return emptyVirus;

}

float VDetect::lambda() const {
    return float(m_currentSize) / float(m_currentCap);
}

float VDetect::deletedRatio() const {
    return float(m_currNumDeleted) / float(m_currentSize);
}

void VDetect::dump() const {

    //cout << "CUR TABLE SIZE: " << getTrueCurrent() << endl;
    //cout << "OLD TABLE SIZE: " << getTrueCurrentOld() << endl;

    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool VDetect::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int VDetect::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) {
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0)
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

int VDetect::getTrueCurrent() const{
    return m_currentSize - m_currNumDeleted;
}

int VDetect::getTrueCurrentOld() const{
    return m_oldSize - m_oldNumDeleted;
}


void VDetect::rehash() {

    //if rehash hasn't been started yet
    if(m_oldTable == nullptr){//creates new currTable and sets previous table to oldTable
        m_oldCap = m_currentCap;//setting variables of old table to curr table in anticipation
        m_oldProbing = m_currProbing;   //of the creation of a new table
        m_oldSize = m_currentSize;
        m_oldNumDeleted = m_currNumDeleted;

        if(m_newPolicy != NONE){
            m_currProbing = m_newPolicy;
        }
        m_currentCap = findNextPrime(m_oldCap * 4);
        Virus *temp = m_currentTable; // temp pointer to store address of m_currentTable
        m_currentTable = new Virus[m_currentCap];//create a new table
        m_oldTable = temp;//point m_oldtable to previous curr-table
        m_currentSize = 0;//resetting variables.
        m_currNumDeleted = 0;
    }

    int bucketsRehashed = (m_oldSize / 4) ; //25% of the nodes
    int transferCounter = 0;//used to count how many nodes were transfered.
    int index = 0;

    while(transferCounter <= bucketsRehashed && index < m_currentCap) {

        if(!m_oldTable[index].m_key.empty() && !(m_oldTable[index] == DELETED)){//if space is filled and it's not deleted
            altInsert(m_oldTable[index]);
            m_oldTable[index] = DELETED;
            m_oldNumDeleted++;
            transferCounter++;
        }
        index++;

        if(getTrueCurrentOld() == 0){//used to exit while-loop if all data has been transferred!
            transferCounter = bucketsRehashed + 1;//there is probably a simpler way to do this, but I confused myself.
        }
    }

    if(getTrueCurrentOld() == 0){//if all the data has been transfered
        delete [] m_oldTable;
        m_oldTable = nullptr;
    }

}

void VDetect::altInsert(Virus virus) {
    //same functionality as insert() but without a reshash() call.
    //used during rehash, inorder to insert nodes from old table to curr table.

    int hashValue = m_hash(virus.getKey()) % m_currentCap;
    int counter = 0;

    //if space is not empty
    //if IDs don't match
    //if space is not deleted
    while(!m_currentTable[hashValue].getKey().empty() && m_currentTable[hashValue].getID() != virus.getID()
          && !(m_currentTable[hashValue] == DELETED)){
        if (m_currProbing == QUADRATIC) {
            hashValue = ((hashValue % m_currentCap) + (counter * counter)) % m_currentCap;
            counter++;
        } else if (m_currProbing == DOUBLEHASH) {
            hashValue = ((hashValue % m_currentCap) + counter * (11 - (hashValue % 11))) % m_currentCap;
            counter++;
        }
    }
    if(m_currentTable[hashValue].getKey().empty() || m_currentTable[hashValue] == DELETED){
        m_currentTable[hashValue] = virus;
        m_currentSize++;
    }

}

ostream& operator<<(ostream& sout, const Virus &virus ) {
    if (!virus.m_key.empty())
        sout << virus.m_key << " (ID " << virus.m_id << ")";
    else
        sout << "";
  return sout;
}

bool operator==(const Virus& lhs, const Virus& rhs){
    return ((lhs.m_key == rhs.m_key) && (lhs.m_id == rhs.m_id));
}