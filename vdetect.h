// CMSC 341 - Spring 2023 - Project 4
#ifndef VDETECT_H
#define VDETECT_H
#include <iostream>
#include <string>
#include "math.h"
using namespace std;
class Grader;   // forward declaration
class Tester;   // forward declaration, will be used for testing
class Virus;    // forward declaration
class VDetect;  // forward declaration
const int MINID = 1000;
const int MAXID = 9999;
const int MINPRIME = 101;   // Min size for hash table
const int MAXPRIME = 99991; // Max size for hash table
#define EMPTY Virus("",0)
#define DELETED Virus("DELETED")
#define DELETEDKEY "DELETED"

typedef unsigned int (*hash_fn)(string);    // declaration of hash function
enum prob_t {NONE, QUADRATIC, DOUBLEHASH};  // types of collision handling policy
#define DEFPOLCY QUADRATIC

const int MAX = 4;
const char ALPHA[MAX] = {'A', 'C', 'G', 'T'};

class Virus{
    public:
    friend class Grader;
    friend class Tester;
    friend class VDetect;
    Virus(string key="", int id=0){m_key = key; m_id = id;}
    string getKey() const {return m_key;}
    int getID() const {return m_id;}
    void setKey(string key){m_key=key;}
    void setID(int id){m_id=id;}
    const Virus& operator=(const Virus& rhs){
        if (this != &rhs){
            m_key = rhs.m_key;
            m_id = rhs.m_id;
        }
        return *this;
    }
    // Overloaded insertion operator
    friend ostream& operator<<(ostream& sout, const Virus &virus );
    // Overloaded equality operator
    friend bool operator==(const Virus& lhs, const Virus& rhs);
    private:
    string m_key;   // the search string used as key in the hash table
    int m_id;       // a unique ID number identifying the object
};

class VDetect{
    public:
    friend class Grader;
    friend class Tester;
    VDetect(int size, hash_fn hash, prob_t probing);
    ~VDetect();
    // Returns Load factor of the new table
    float lambda() const;
    // Returns the ratio of deleted slots in the new table
    float deletedRatio() const;
    // insert only happens in the new table
    bool insert(Virus virus);
    // remove can happen from either table
    bool remove(Virus virus);
    // find can happen in either table
    Virus getVirus(string key, int id) const;
    // request a change in collision handling policy
    void changeProbPolicy(prob_t policy);
    // dumps the contents of the two tables
    void dump() const;

    private:
    hash_fn    m_hash;          // hash function
    prob_t     m_newPolicy;     // stores the change of policy request

    Virus*     m_currentTable;  // hash table
    int        m_currentCap;    // hash table size (capacity)
    int        m_currentSize;   // current number of entries
                                // m_currentSize includes deleted entries 
    int        m_currNumDeleted;// number of deleted entries
    prob_t     m_currProbing;   // collision handling policy

    Virus*     m_oldTable;      // hash table
    int        m_oldCap;        // hash table size (capacity)
    int        m_oldSize;       // current number of entries
                                // m_oldSize includes deleted entries
    int        m_oldNumDeleted; // number of deleted entries
    prob_t     m_oldProbing;    // collision handling policy

    //private helper functions
    bool isPrime(int number);
    int findNextPrime(int current);

    /******************************************
    * Private function declarations go here! *
    ******************************************/

    int getTrueCurrent() const;
    int getTrueCurrentOld() const;
    void rehash();//does the rehash operation
    void altInsert(Virus virus);//does the same thing as insert, but doesn't check for rehash.
                                    //used in rehash function to insert virus from old to curr.

};
#endif