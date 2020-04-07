#include "addgraphform.h"
#include "ui_addgraphform.h"

#include "mainwindow.h"
#include <QMessageBox>

#include <iostream>

AddGraphForm::AddGraphForm(QWidget *parent, int _index, int tabNumber, int _CreatedGraphsCounter, QColor _color) :
    QWidget(parent),
    ui(new Ui::AddGraphForm),
    index(_index),
    TabNumber(tabNumber),
    color(_color),
    CreatedGraphsCounter(_CreatedGraphsCounter)
{
    ui->setupUi(this);
    ui->label->setText(QString("Formula %0").arg(CreatedGraphsCounter));

    ui->ColorLabel->setAutoFillBackground(true);
    QPalette palette = ui->ColorLabel->palette();
    palette.setColor(ui->ColorLabel->backgroundRole(), color);
    ui->ColorLabel->setPalette(palette);

    std::cout<<"AddGraphForm(): _index="<<_index<<" index="<<index
             <<" index+1="<<index+1<<std::endl;

    //I DO NOT WHY THIS WORKS
    //IF TO UNCOMMNET IT, WHEN PRESS ON NO BUTTON IN on_DeleteButton_clicked(),
    // the QMessageBox::warning(...) message is called twice.

    //connect(ui->DeleteButton, SIGNAL(clicked()), this, SLOT(on_DeleteButton_clicked()));

}

AddGraphForm::~AddGraphForm()
{
    delete ui;
}

QString AddGraphForm::GetFormula()
{
    return ui->lineEdit->text();
}



void AddGraphForm::on_DeleteButton_clicked()
{
    QMessageBox::StandardButton reply =
            QMessageBox::warning(this, "Delete Graph",
                                 QString("Delete Graph %0").arg(CreatedGraphsCounter),
                                 QMessageBox::Yes | QMessageBox::No);
    if(reply == QMessageBox::Yes) emit DeleteRequested();
    else if(reply == QMessageBox::No){return;}
}
