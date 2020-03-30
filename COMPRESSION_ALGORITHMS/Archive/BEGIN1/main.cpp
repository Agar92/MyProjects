#include <queue>
#include <iostream>
#include <sstream>
using namespace std;
 
struct Node {
    string value, code;
    unsigned amount;
    Node * left;
    Node * right;
    //компаратор
    bool operator() (const Node& x, const Node& y) const {
        return x.amount > y.amount;
    }
    //конструктор по умолчанию нужен для создания объекта-компаратора
    Node (const string& value = "", unsigned amount = 0, Node * left = 0, Node * right = 0) {
        this->value = value;   //множество символом узла
        this->code = "";       //строковое представление битового кода узла
        this->amount = amount; //сколько раз втретилось
        this->left = left;     //левый ребенок
        this->right = right;   //правый ребенок
    }
    //генерация DOT-описания полученного дерева 
    string to_str() {
        ostringstream x;
        if (left != 0 || right != 0) { // дерево таково, что либо есть оба ребенка либо нет ни одного
            x << "\t\"" << code << ": " << value << "[" << amount << "]\" -> \""  << left->code << ": " << left->value << "[" << left->amount << "]\";\n";
            x << left->to_str();
            x << "\t\"" << code << ": " << value << "[" << amount << "]\" -> \""  << right->code << ": " << right->value << "[" << right->amount << "]\";\n";
            x << right->to_str();
        } else {
            x << "\t\"" << code << ": " << value << "[" << amount << "]\" [shape=box, style=filled, fillcolor=green];\n";
        }
        return x.str();
    }
    //  объединение деревьев
    Node * join (Node x) {
        return new Node(x.value + value, x.amount + amount, new Node(x), this);
    }
    // проход по дереву с генерацией кода
    void traversal_code(string code) {
        this->code = code;
        if (left != 0 || right != 0) {
            left->traversal_code(code + "1");
            right->traversal_code(code + "0");
        }
    }
    // строим дерево по алгоритму Хаффмана
    static Node * builder(priority_queue<Node, vector<Node>, Node> graph) {
        while (graph.size() > 1) {
            Node *n = new Node(graph.top());
            graph.pop();
            graph.push(*n->join(*new Node(graph.top())));
            graph.pop();
        }
        return new Node(graph.top());
    }
};
 
unsigned amounts[256]; // массив счетчиков встречаемости символов
int main() {
    string s;
    getline (std::cin, s); // читаем строку вместе с пробелами
    for(auto i: s) amounts[i]++;
    priority_queue<Node, vector<Node>, Node> graph;
    for(int i = 'a'; i <= 'z'; i++) // записываем в очередь с приоритетами 
        if (amounts[i] > 0) graph.emplace(s=(char)i, amounts[i]);
    Node *n = Node::builder(graph);
    n->traversal_code("");
    // генерируем ссылку на сервис Google для генерации изображений по описанию графа
    cout << "http://chart.apis.google.com/chart?cht=gv&chl=" << endl;
    // генерируем DOT-описание полученного дерева для отрисовки
    cout << "Digraph G {\n" << n->to_str() << "}";
    // Если вывод программы скопировать и вставить в адресную строку браузера, то увидим картинку.  
}
