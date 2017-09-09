//
// This file contains implementation of  several functions
// to test the Set class.
//

#include "Tests.h"
#include "Set.h"

extern double calculateTVD(
    vector<pair<int, double>> rgA,
    vector<pair<int, double>> rgB
    );

void TestLengthDistribution()
{
    char *rgStrA[] = {"a", "aab", "caab", "abababa", "aba", NULL};
    char *rgStrB[] = {"ab", "aab", "daa", "abbaba", "cba", NULL};

    vector<pair<int, double>> rgA, rgB;
    Set A, B, C;
    int i;

    for (i = 0; rgStrA[i]; i++)
    {
        A.AddElement(rgStrA[i]);
    }
    A.Print();
    rgA = A.LengthDistribution();
    for (i = 0; i < (int)rgA.size(); i++)
    {
        printf("%d: %d\n", rgA[i].first, rgA[i].second);
    }

    for (i = 0; rgStrB[i]; i++)
    {
        B.AddElement(rgStrB[i]);
    }
    B.Print();
    rgB = B.LengthDistribution();
    for (i = 0; i < (int)rgB.size(); i++)
    {
        printf("%d: %d\n", rgB[i].first, rgB[i].second);
    }

    printf("TVD = %0.5lf\n", calculateTVD(rgA, rgB));

}

void TestSymmetricDifference()
{
    char *rgStrA[] = {"aaa", "aab", "caa", "abb", "aba", NULL};
    char *rgStrB[] = {"aaa", "aab", "daa", "abb", "cba", NULL};
    //char *rgStrB[] = {NULL};

    Set A, B, C;
    int i;

    for (i = 0; rgStrA[i]; i++)
    {
        A.AddElement(rgStrA[i]);
    }
    A.Print();

    for (i = 0; rgStrB[i]; i++)
    {
        B.AddElement(rgStrB[i]);
    }
    B.Print();

    C = A.SymmetricDifference(B);
    C.Print();
}


void TestDifference()
{
    char *rgStrA[] = {"aaa", "aab", "caa", "abb", "aba", NULL};
    char *rgStrB[] = {"aaa", "aab", "daa", "abb", "cba", NULL};
    //char *rgStrB[] = {NULL};

    Set A, B, C;
    int i;

    for (i = 0; rgStrA[i]; i++)
    {
        A.AddElement(rgStrA[i]);
    }
    A.Print();

    for (i = 0; rgStrB[i]; i++)
    {
        B.AddElement(rgStrB[i]);
    }
    B.Print();

    C = A.Difference(B);
    C.Print();
}


void TestIntersection()
{
    char *rgStrA[] = {"aaa", "aab", "caa", "abb", "aba", NULL};
    //char *rgStrB[] = {"aaa", "aab", "daa", "abb", "cba", NULL};
    char *rgStrB[] = {NULL};

    Set A, B, C;
    int i;

    for (i = 0; rgStrA[i]; i++)
    {
        A.AddElement(rgStrA[i]);
    }
    A.Print();

    for (i = 0; rgStrB[i]; i++)
    {
        B.AddElement(rgStrB[i]);
    }
    B.Print();

    C = A.Intersection(B);
    C.Print();
}


void TestUnion()
{
    char *rgStrA[] = {"aaa", "aab", "caa", "abb", "aba", NULL};
    char *rgStrB[] = {"aaa", "aab", "daa", "abb", "cba", NULL};

    Set A, B, C;
    int i;


    A.StartBulkAddElements();
    for (i = 0; rgStrA[i]; i++)
    {
        A.AddElement(rgStrA[i]);
    }
    A.EndBulkAddElements();
    A.Print();

    B.StartBulkAddElements();
    for (i = 0; rgStrB[i]; i++)
    {
        B.AddElement(rgStrB[i]);
    }
    B.EndBulkAddElements();
    B.Print();

    C = A.Union(B);
    C.Print();
}

void TestAddElement()
{
    Set A;

    A.AddElement("aaa"); // aaa
    A.Print();

    A.AddElement("aab"); // aaa, aab
    A.Print();

    A.AddElement("caa"); // aaa, aab, caa
    A.Print();

    A.AddElement("abb"); // aaa, aab, abb, caa
    A.Print();

    A.AddElement("aba"); // aaa, aab, aba, abb, caa
    A.Print();
}
