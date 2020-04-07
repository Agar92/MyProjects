#ifndef ADDGRAPHFORM_H
#define ADDGRAPHFORM_H

#include <QWidget>

namespace Ui {
class AddGraphForm;
}

class AddGraphForm : public QWidget
{
    Q_OBJECT

public:
    explicit AddGraphForm(QWidget *parent = nullptr, int _index=0, int tabNumber=0, int _CreatedGraphsCounter=0, QColor _color = Qt::black);
    ~AddGraphForm();
    int GetIndex(){return index;}
    int GetCountNumber(){return CreatedGraphsCounter;}
    int GetTabNumber(){return TabNumber;}
    void SetTabNumber(int _tabNumber){TabNumber=_tabNumber;}
    QString GetFormula();
    QColor GetColor(){return color;}

private slots:
    void on_DeleteButton_clicked();

signals:
    void DeleteRequested();

private:
    Ui::AddGraphForm *ui;
    int index;
    int TabNumber;
    QColor color;
    int CreatedGraphsCounter;
};

#endif // ADDGRAPHFORM_H
