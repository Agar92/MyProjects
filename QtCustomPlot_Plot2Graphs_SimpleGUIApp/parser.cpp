#include "parser.h"

set<char> DelimSet;
map <string, int> prior;
Ops ops;

//Функция считывает выражение в строку "expr" и ищем переменные или массивы
void ReadExpressionFromStream(ifstream &inp, string &expr, Variables &var, Massives &mas) {
    getline(inp, expr);
    string temp;
    int pos;
    while (!inp.eof())
    {
        getline(inp, temp);
        //Если встретили '=', то это переменная, заносим ее имя и значение в массив
        pos = temp.find('=');
        if(pos>0)
        {
            string name = temp.substr(0, pos);
            double value = atof(temp.substr(pos + 1).c_str());
            var[name] = value;
        }
        //Если нашли '{', это массив, заносим имя массива и значения в соответствующие массивы
        pos = temp.find('{');
        if(pos>0)
        {
            string name = temp.substr(0,pos);
            //Ищем значения массива и запоминаем их
            int pos1=pos, pos2;
            do {
                pos2 = temp.find(',');
                double value = atof(temp.substr(pos1+1, pos2).c_str());
                mas[name].push_back(value);
                if(pos2 == -1) break;
                temp[pos2] = ' ';
                pos1 = pos2;
            }while (pos2 > 0);
        }
    }
    return;
}


//Разбиваем выражение на токены
void CreateTokensFromExpression(string &expr, tokens &texpr) {
    CreateSetOfDelimiters();
    string ex = expr + " ";
    string name;

    //Получаем имя токена
    int i = 0;
    while (i < ex.size() - 1)
    {
        name = "";
        //Если текущий символ разделитель
        if (IsDelimiter(ex[i]))
        {
            if (ex[i] == ' ')
            { //Пробел просто перепрыгиваем
                i++;
                continue;
            }
            name = ex[i]; //Любой другой добавляем в имя токена
            i++;
        }
        else
        {
            while (!IsDelimiter(ex[i]))
            {
            /*Если не разделитель, например, переменная или имя массива,
            Считываем его полностью */
                name += ex[i];
                i++;
            }
        }
        //Заносим получившийся токен в список токенов
        texpr.push_back(token(name, var));
    }

    //Раздаем получившимся токенам типы
    for (int j = 0; j < texpr.size(); j++)
    {
        if (texpr[j].name[0] == '[') {
            texpr[j].type = op_sbr;
            continue;
        }
        if (texpr[j].name[0] == ']')
        {
            texpr[j].type = cl_sbr;
            continue;
        }
        if (texpr[j].name[0] == '(')
        {
            texpr[j].type = op_br;
            continue;
        }
        if (texpr[j].name[0] == ')')
        {
            texpr[j].type = cl_br;
            continue;
        }
        if (isdigit(texpr[j].name[0])) {
            texpr[j].type = num;
            continue;
        }

        //mas
        if (isalpha(texpr[j].name[0])) {
            if(j < texpr.size() - 1 && texpr[j+1].name[0] == '[')
                texpr[j].type = mas;
            //continue;
        }

        if (isalpha(texpr[j].name[0])) {
            if (j < texpr.size() - 1 && texpr[j+1].name[0] == '(')
                texpr[j].type = func;
            continue;
        }

        texpr[j].type = op;
    }
    //Проверяем минус и !, что это префиксные операции
    for (int j = 0; j < texpr.size(); j++)
    {
        if (texpr[j].name == "-" && (j == 0 || texpr[j - 1].type == op_br))
            texpr[j].name = "opposite";
        if (texpr[j].name == "!" && (j == texpr.size()-1 || texpr[j+1].type == cl_br || texpr[j+1].type == op))
            texpr[j].name = "factorial";
    }
    return;
}

//Переводим выражение в постфиксную запись
void CreatePostfixFromTokens(tokens &texpr, tokens &pexpr) {
    //Задаем приоритеты операций
    CreatePrior();
    stack <token> TStack;

    //Ловим токены и работаем по алгоритму
    for (int i = 0; i < texpr.size(); i++) {
        switch (texpr[i].type) {
            case var:
            case num:
                pexpr.push_back(texpr[i]);
                break;

            case op_br:
                TStack.push(texpr[i]);
                break;

            case cl_br:
                while (TStack.top().type != op_br) {
                    pexpr.push_back(TStack.top());
                    TStack.pop();
                }
                TStack.pop();
                break;

            case op_sbr:
                TStack.push(texpr[i]);
                break;

            case cl_sbr:
                while (TStack.top().type != op_sbr) {
                    pexpr.push_back(TStack.top());
                    TStack.pop();
                }
                TStack.pop();
                break;

            case op:
                if (TStack.size()) {
                    while (TStack.size() && ((TStack.top().type == op && prior[texpr[i].name] <= prior[TStack.top().name]) ||
                        TStack.top().type == func ||
                        TStack.top().type == mas)) {
                        pexpr.push_back(TStack.top());
                        TStack.pop();
                    }
                }
                TStack.push(texpr[i]);
                break;

            case mas:
                while (TStack.size() && TStack.top().type == mas) {
                    pexpr.push_back(TStack.top());
                    TStack.pop();
                }
                TStack.push(texpr[i]);
                break;

            case func:
                while (TStack.size() && TStack.top().type == func) {
                    pexpr.push_back(TStack.top());
                    TStack.pop();
                }
                TStack.push(texpr[i]);
                break;
        }
    }

    while (TStack.size()) {
        pexpr.push_back(TStack.top());
        TStack.pop();
    }

    return;
}



//Считаем результат выражения
double ResultExpr(tokens &pexpr, Variables &expvars, Massives &varmas) {
    CreateOps();
    stack <double> s;

    for (int i=0; i<pexpr.size(); i++) {
        switch(pexpr[i].type) {
            case num: {
                //stdlib.h functions somewhy do not work
                //the return a=0 instead of a=0.5
                //and a=1 instead of a=1.5
                //s.push(atof(pexpr[i].name.c_str()));
                //That is why use this:
                string STRING(pexpr[i].name.c_str());
                istringstream str(STRING);
                double val;
                str>>val;
                s.push(val);

                /*
                cout<<"pexpr[i].name.c_str()="
                    <<pexpr[i].name.c_str()
                    <<" atof(pexpr[i].name.c_str())="
                    <<atof(pexpr[i].name.c_str())
                    <<endl;
                */
                }
                break;

            case var: {
                    Variables::iterator Vit;
                    for(Vit=expvars.begin(); Vit!=expvars.end(); Vit++) {
                        if(Vit->first == pexpr[i].name) {
                            s.push(Vit->second);
                            break;
                        }
                    }
                }
                break;

            case func:
            case op: {
                    Ops::iterator Oit;
                    for(Oit=ops.begin(); Oit!=ops.end(); Oit++) {
                        if(Oit->first == pexpr[i].name) {
                            s.push(Oit->second(s));
                            break;
                        }
                    }
                }
                break;

            case mas: {
                int index = s.top();
                s.pop();

                Massives::iterator it;
                for(it = varmas.begin(); it != varmas.end(); it++) {
                    if(it->first == pexpr[i].name)
                        s.push(it->second[index]);
                }
            }
        }
    }

    return s.top();
}


