#include <iostream>
using namespace std;
int main()
{
    int a[10], b[10];
    int mn, index, i, j, temp, m, n;
    a[0] = 1;
    a[1] = 7;
    a[2] = 3;
    a[3] = -1;
    a[4] = 6;
        for ( i = 0 ; i < 5 ; i++ )
        {
            mn = a[i];
            index = i;
            for ( j = 0 ; j < 5 ; j++ )
            {
                if ( mn > a[j] && a[j] != 999999999 )
                {
                    mn = a[j];
                    index = j;
                }
            }
            a[index] = 999999999;
            b[i] = index;
            cout << "Index : " << index << endl;
            for ( m = 0 ; m < 5 ; m++ )
                cout << a[m] << " ";
                cout << '\n';
            for ( n = 0 ; n < 5 ; n++ )
                cout << b[n] << " ";
            cout << '\n';
        }


}
