// This is an over-engineered example of a Fisher-Yates algorithm
// but it should work quite robustly also generates the list of 
// equal numbers of each quantum number. 

// Importantly template classes cannot return arrays and so we use a class 
// structure instead of a template function which can only return the typename 
// template object 


#ifndef __RANDOM_STATE
#define __RANDOM_STATE

#include <iostream> 
#include <random> 
#include <vector> 
#include <iterator>
#include <array>
#include <string>
#include <tuple>
#include <stdexcept>

template<typename Iter, typename RandomGenerator> 
Iter select_randomly(Iter start, Iter end, RandomGenerator& g){
    std::uniform_int_distribution<> dis(0,std::distance(start,end)-1);
    std::advance(start,dis(g)); 
    return start;
}

//Randomly returns a pointer (iterator) to a variable between the start and end of the list
//Using the template class enables us to swap between different variable types 
template<typename Iter> 
Iter select_randomly(Iter start, Iter end) { 
    static std::random_device rd; 
    static std::mt19937 gen(rd()); 
    return select_randomly(start,end,gen); 
}

//Takes in a std::string {"Up","Down"} and randomly selects either one of them
//Size gives the number of elements in the array from which we are randomly selecting 
//N gives the number of selections required which in the case for runauto.cc  is the 
//number of states in the unit cell

//Symmetric contributions only and therefore we need to have equal numbers of each 
//quantum number contribution and therefore N%len == 0 otherwise throw an error 

template <typename T, int size> 
class randomStatesSz{

    public: 
    
        int len;
        int i = 0;
        std::vector<int> numberList;
        std::vector<T> finalString;
        std::vector<std::tuple<T,int>> counterArray;


    randomStatesSz(T (&array)[size], int N){

        T* initString[N]; 

        len = sizeof(array)/sizeof(array[0]);

        if ((N%len) != 0){
            throw std::invalid_argument("N%lenArray == 0");
        }
        
        while (i != len){

            for (int j = i*floor(N/len); j < (i+1)*floor(N/len); j++){
                initString[j] = &array[i]; 
            }
            i += 1; 
            
        }

        for (int i = 0; i < N; i++){
            numberList.push_back(i); 
        }

        while (numberList.size() > 0){
            auto selected = select_randomly(numberList.begin(), numberList.end()); 
            int intSelected = *selected;
            finalString.push_back(*initString[intSelected]);
            numberList.erase(selected);
            i++;
        }

    }

    
    std::vector<T> getRandomNumbers() {return finalString;};

    std::vector<std::tuple<T,int>> randomTotals(){

        for (T i : finalString){

            int counter = 0; 

            for (std::tuple<T,int> &j : counterArray){

                if (i==std::get<0>(j)){
                    std::get<1>(j) += 1;
                    counter +=1; 
                }

            }

            if (counter == 0){
                counterArray.push_back(std::make_tuple(i,1));
            }
            
        }

        return counterArray;

    }
         

         
};


/// Absolutely no reason to do this by pointers bar the fact that i like them 
/// This is equivalent to just randomly selecting a number, removing it and the 
/// element in the number. So a selection with replacement. Would be easier but 
/// at least this way we have something which is horribly complicated and over
/// engineered

#endif