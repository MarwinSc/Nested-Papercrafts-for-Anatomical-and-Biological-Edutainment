#include <stdio.h>
#include <iostream>
using namespace std;

#ifdef __cplusplus
extern "C" int hello(int i)
#else
int hello(int i)
#endif
{
    cout << "Yo";
    return i*2;
}