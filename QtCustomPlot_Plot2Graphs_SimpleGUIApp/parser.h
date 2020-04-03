#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <map>
#include <set>
#include <vector>
#include <stack>
#include <iostream>
#include <fstream>
#include <cmath>

#include <stdlib.h>
#include <cstdlib>
#include <sstream>

using namespace std;

//Массив переменных
typedef map<string, double> Variables;
//Массив "массивов"
typedef map<string, vector<double> > Massives;

//Функция считывает выражение в строку "expr" и ищем переменные или массивы
void ReadExpressionFromStream(ifstream &inp, string &expr, Variables &var, Massives &mas);

// Типы токенов
enum tokentype {
    //Переменная, константа, (, ), функция, операция, массив, {, }
    var, num, op_br, cl_br, func, op, mas, op_sbr, cl_sbr
};
// Структура токена
struct token {
    string name;
    tokentype type;

    //Конструкторы
    token(string str, tokentype typ) {
        name = str;
        type = typ;
    }
    token() {}
};
//Массив токенов
typedef vector<token> tokens;
//Множество разделителей
extern set<char> DelimSet;
//Разделители
const string delimiters = " ()+/*-^&|!%[]";
//Инициализирует множество разделителей
inline void CreateSetOfDelimiters()
{
    for (int i = 0; i < delimiters.size(); i++)
        DelimSet.insert(delimiters[i]);
    return;
}

//Проверка, является ли символ разделителем
inline bool IsDelimiter(char sym)
{
    return DelimSet.count(sym) > 0;
}

//Разбиваем выражение на токены
void CreateTokensFromExpression(string &expr, tokens &texpr);


//Приоритеты операций
extern map <string, int> prior;
//Функция выставляет приоритеты операций
inline void CreatePrior()
{
    prior["+"] = 10;
    prior["-"] = 10;
    prior["*"] = 20;
    prior["/"] = 20;
    prior["^"] = 30;
    prior["opposite"] = 10;
    prior["factorial"] = 30;
    prior["%"] = 20;
    prior["&"] = 5;
    prior["|"] = 5;
    prior["!"] = 40;
}

//Переводим выражение в постфиксную запись
void CreatePostfixFromTokens(tokens &texpr, tokens &pexpr);


//Указатель на функцию(для операций)
typedef double(*func_type)(stack<double> &);

//Массив операций
typedef map<string, func_type> Ops;
extern Ops ops;



//Реализация доступных операций
inline double fact(double n)
{
    if(n == 0)
        return 1;
    return n*fact(n-1);
}
inline double op_plus(stack <double> &s)
{
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return a+b;
}
inline double op_minus(stack <double> &s)
{
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return b-a;
}
inline double op_mul(stack <double> &s)
{
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return a*b;
}
inline double op_div(stack <double> &s)
{
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return b/a;
}
inline double op_deg(stack <double> &s)
{
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    //b^a!!
    cout<<"a="<<a<<" b="<<b<<" pow(b,a)="
        <<pow(b,a)<<endl;
    return pow(b,a);
}
inline double op_opposite(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return -a;
}
inline double op_factorial(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return fact(a);
}
inline double op_odiv(stack <double> &s)
{
    long long a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return b%a;
}
inline double op_and(stack <double> &s) {
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return a&&b;
}
inline double op_or(stack <double> &s)
{
    double a,b;
    a = s.top();
    s.pop();
    b = s.top();
    s.pop();
    return a||b;
}
inline double op_not(stack <double> &s) {
    double a;
    a = s.top();
    s.pop();
    return !a;
}
inline double op_sin(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return sin(a);
}
inline double op_cos(stack <double> &s) {
    double a;
    a = s.top();
    s.pop();
    return cos(a);
}
inline double op_tan(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return tan(a);
}
inline double op_asin(stack <double> &s) {
    double a;
    a = s.top();
    s.pop();
    return asin(a);
}
inline double op_acos(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return acos(a);
}
inline double op_atan(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return atan(a);
}
inline double op_exp(stack <double> &s) {
    double a;
    a = s.top();
    s.pop();
    return exp(a);
}
inline double op_log(stack <double> &s)
{
    double a;
    a = s.top();
    s.pop();
    return log(a);
}

//Инициализация массива операций
inline void CreateOps()
{
    ops["+"] = op_plus;
    ops["-"] = op_minus;
    ops["*"] = op_mul;
    ops["/"] = op_div;
    ops["^"] = op_deg;
    ops["opposite"] = op_opposite;
    ops["factorial"] = op_factorial;
    ops["%"] = op_odiv;
    ops["&"] = op_and;
    ops["|"] = op_or;
    ops["!"] = op_not;
    ops["sin"] = op_sin;
    ops["cos"] = op_cos;
    ops["tan"] = op_tan;
    ops["acos"] = op_acos;
    ops["asin"] = op_asin;
    ops["atan"] = op_atan;
    ops["exp"] = op_exp;
    ops["log"] = op_log;
    return;
}

//Считаем результат выражения
double ResultExpr(tokens &pexpr, Variables &expvars, Massives &varmas);




#endif // PARSER_H
