#include <QCoreApplication>
#include <iostream>

#include "Globals.h"
#include "parser.h"

extern const int N;

using namespace std;

//From http://mindhalls.ru/symb-expr/.

int main() {

    string expression = "x^6*log(x)\nx=3";
    ofstream out("test.txt", ios::out);
    out<<expression;
    out.close();


    tokens texpr, pexpr;
    Variables expvars;
    Massives expmasvars;
    string expr;
    ifstream file("test.txt");

    ReadExpressionFromStream(file, expr, expvars, expmasvars);

    cout << "Expr:" << endl;
    cout << expr << endl;

    Variables::iterator it;
    for (it = expvars.begin(); it != expvars.end(); it++)
        cout << it->first << '=' << it->second << endl;

    Massives::iterator it1;
    for (it1 = expmasvars.begin(); it1 != expmasvars.end(); it1++) {
        cout << it1->first << '{';
        for(int i=0; i<it1->second.size(); i++) {
            if(i == it1->second.size()-1)
                cout << it1->second[i];
            else
                cout << it1->second[i] << ',';
        }
        cout << '}' << endl;
    }
    cout << endl;


    CreateTokensFromExpression(expr, texpr);

    cout << "Token:" << endl;
    for (int i = 0; i < texpr.size(); i++)
        cout << texpr[i].name << ' ';
    cout << endl << endl;

    CreatePostfixFromTokens(texpr, pexpr);

    cout << "Pexpr:" << endl;
    for (int i = 0; i < pexpr.size(); i++)
        cout << pexpr[i].name << ' ';
    cout << endl << endl;

    cout << "Result:" << endl;
    cout << ResultExpr(pexpr, expvars, expmasvars) << endl;

    return 0;
}
