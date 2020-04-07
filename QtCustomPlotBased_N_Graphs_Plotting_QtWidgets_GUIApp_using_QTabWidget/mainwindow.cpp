#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "qcustomplot.h"
#include "addgraphform.h"

#include <iostream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //Скрываем элементы управления и рамку окна...
    this->setWindowFlags(Qt::FramelessWindowHint);
    //Включаем прозрачность главной формы...
    setAttribute(Qt::WA_TranslucentBackground );

    AddGraphTabsNumber=0;
    ColorCounter=0;
    CreatedGraphsCounter=0;
    for(int i=0; i<NUMBER_OF_PLOTS; ++i) alive[i]=false;
//    connect(form, SIGNAL(on_DeleteButton_clicked()), this, SLOT(on_tabWidget_Close()));

    //ui->tabWidget->setCurrentIndex(1);
    settingsForm = new SettingsForm(nullptr);
    connect(settingsForm, SIGNAL(SettingsChanged(QVector<double> &)), this, SLOT(ApplySettingChanges(QVector<double> &)));
    ui->tabWidget->addTab(settingsForm, QString("Settings"));

    ui->tabWidget->setCurrentIndex(0);

    Xmin = 0.1;
    Xmax = 100.0;
    Ymin =-100.0;
    Ymax = 100.0;
    ui->widget->xAxis->setRange(Xmin, Xmax);//Для оси Ox
    ui->widget->yAxis->setRange(Ymin, Ymax);//Для оси Oy
    //И перерисуем график на нашем widget
    ui->widget->replot();
}

MainWindow::~MainWindow()
{
    delete ui;
    delete settingsForm;
}

void MainWindow::PlotGraphs()
{
    QVector<int> CountNumber;
    QVector<QString> FormulaVector;
    QVector<QColor> ColorVector;
    QVector<int> IndexVector;
    for(int i=0; i<NUMBER_OF_PLOTS; ++i)
    {
        if(alive[i])
        {
            qDebug() << form[i]->GetFormula()
                     << "   " << form[i]->GetFormula().size()
                     << endl;

            if(form[i]->GetFormula().size()>0)
            {
                CountNumber.push_back(form[i]->GetCountNumber());
                FormulaVector.push_back(form[i]->GetFormula());
                ColorVector.push_back(form[i]->GetColor());
                IndexVector.push_back(form[i]->GetIndex());
            }

        }
    }

    //number of points plotted:
    const int N=10000;
    QVector<double> x(N), y(N);
    //fill Ox:
    for(int i=0; i<N; ++i) x[i]=Xmin+i * (Xmax-Xmin)/(N-1);
    ui->widget->clearGraphs();//Если нужно, но очищаем все графики

    //prepare y_i values:
    for(int I=0; I<FormulaVector.size(); ++I)
    {
        tokens texpr, pexpr;
        Variables expvars;
        Massives expmasvars;
        string expr = FormulaVector.at(I).toUtf8().constData();

        std::cout<<"I="<<I<<" expr="<<expr<<" size="<<expr.size()
                 <<std::endl;

        if(expr.size()==0)
        {
            QMessageBox::warning(this, "Incorret input", "Empty input");
            /*
            QMessageBox msgBox;
            msgBox.setText("Empty input!");
            //msgBox.setIcon(QMessageBox::Information);
            //msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            //QTimer::singleShot(1000, &msgBox, SLOT(close()));
            //msgBox.close();
            */
        }
        else
        {
            CreateTokensFromExpression(expr, texpr);
            CreatePostfixFromTokens(texpr, pexpr);
            for(int i=0; i<N; ++i)
            {
                expvars["x"]=x[i];
//                cout << "Result:" << endl;
                y[i]=ResultExpr(pexpr, expvars, expmasvars);
//                cout << y.at(i) << endl;
            }

            //Добавляем один график в widget
            ui->widget->addGraph();
            //Говорим, что отрисовать нужно график по нашим двум массивам x и y
            ui->widget->graph(I)->setData(x, y);
            QPen pen;  // создаёт перо по умолчанию
            pen.setStyle(Qt::SolidLine);//(Qt::DashDotLine);
            pen.setWidth(3);
            pen.setBrush(form[(IndexVector.at(I))]->GetColor());
            ui->widget->graph(I)->setPen(pen);
            //Подписываем оси Ox и Oy
            ui->widget->xAxis->setLabel("x");
            ui->widget->yAxis->setLabel("y");
            double minX = Xmin, maxX = Xmax;
            //Установим область, которая будет показываться на графике
            ui->widget->xAxis->setRange(minX, maxX);//Для оси Ox
            //Для показа границ по оси Oy сложнее, так как надо правильно
            //вычислить минимальное и максимальное значение в векторах:
            double minY = Ymin;//y[0]
            double maxY = Ymax;//y[0];
            ui->widget->yAxis->setRange(minY, maxY);//Для оси Oy
        }
    }
    //И перерисуем график на нашем widget
    ui->widget->replot();
    //QTimer::singleShot(2000, this, SLOT(ui->progressBar->setValue(0)));
}

void MainWindow::on_tabWidget_Close()
{
    std::cout<<"HERE"<<std::endl;
    int IND = ui->tabWidget->currentIndex() - 2 + 1;

    int INDEX=0;
    int TABNUMBER=0;
    for(int i=0; i<NUMBER_OF_PLOTS; ++i)
    {
        if(alive[i])
            std::cout<<"ui->tabWidget->currentIndex()="<<ui->tabWidget->currentIndex()
                     <<" form[i]->GetTabNumber()="<<form[i]->GetTabNumber()
                     <<" form[i]->GetIndex()="<<form[i]->GetIndex()
                     <<" alive["<<i<<"]="<<alive[i]<<std::endl;
        else
            std::cout<<"ui->tabWidget->currentIndex()="<<ui->tabWidget->currentIndex()
                     <<" form[i]->GetTabNumber()="<<"No"
                     <<" form[i]->GetIndex()="<<"No"
                     <<" alive["<<i<<"]="<<alive[i]<<std::endl;

        if(alive[i] == true && form[i]->GetTabNumber() == ui->tabWidget->currentIndex())
        {
            INDEX = form[i]->GetIndex();
            TABNUMBER = form[i]->GetTabNumber();
        }
    }

    for(int i=0; i<NUMBER_OF_PLOTS; ++i)
    {
        if(i != INDEX && alive[i] == true)
        {
            if(form[i]->GetTabNumber()>=TABNUMBER) form[i]->SetTabNumber(form[i]->GetTabNumber()-1);
        }
    }


    int Index = form[INDEX]->GetIndex();

    std::cout<<"on_tabWidget_Close(): Index="<<Index
             <<" IND="<<IND
             <<" curTab="<<ui->tabWidget->currentIndex()
             <<" Tab count="<<ui->tabWidget->count()
             <<" AddGraphTabsNumber="<<AddGraphTabsNumber
             <<" currentIndex="<<ui->tabWidget->currentIndex()
             <<" INDEX="<<INDEX
             <<std::endl;

    std::cout<<"on_tabWidget_Close(): ui->tabWidget->currentIndex()="<<ui->tabWidget->currentIndex()
             <<std::endl;

    ui->tabWidget->removeTab(ui->tabWidget->currentIndex());
    delete form[INDEX];
    alive[INDEX] = false;
    --AddGraphTabsNumber;
}

void MainWindow::ApplySettingChanges(QVector<double> & vec)
{
    std::cout<<"Before settingsForm->GetSettings():"<<std::endl;
    //QVector<double> vec = settingsForm->GetSettings();

    std::cout<<"SIZE="<<vec.size()<<std::endl;
    std::cout<<"vec: "<<vec[0]<<" "
             <<vec[1]<<" "<<vec[2]
             <<" "<<vec[3]<<std::endl;

    Xmin = vec.at(0);
    Xmax = vec.at(1);
    Ymin = vec.at(2);
    Ymax = vec.at(3);

    std::cout<<"Xmin="<<Xmin<<" Xmax="<<Xmax
             <<" Ymin="<<Ymin<<" Ymax="<<Ymax
             <<std::endl;

    PlotGraphs();

}

void MainWindow::on_AddGraphButton_clicked()
{
    if(AddGraphTabsNumber<NUMBER_OF_PLOTS)
    {
        int index = ui->tabWidget->count() - 2 + 1;
        if(ColorCounter>8) ColorCounter=0;
        QColor color = COLORS_STORED[ColorCounter];
        ++ColorCounter;
        ++CreatedGraphsCounter;

        form[index-1] = new AddGraphForm(nullptr, index-1, AddGraphTabsNumber + 2, CreatedGraphsCounter, color);
        alive[index-1] = true;

        std::cout<<"on_AddGraphButton_clicked(): currentIndex="<<AddGraphTabsNumber + 2
                 <<std::endl;

        ++AddGraphTabsNumber;

        std::cout<<"QQQ: index="<<index<<" index-1="<<index-1
                 <<" AddGraphTabsNumber="<<AddGraphTabsNumber
                 <<std::endl;

//        connect(form, SIGNAL(form->DeleteButton->clicked()), this, SLOT(this->on_tabWidget_Close()));

        connect(form[index-1] , SIGNAL(DeleteRequested()), this, SLOT(on_tabWidget_Close()));

        ui->tabWidget->addTab(form[index-1], QString("Graph %0").arg(CreatedGraphsCounter));
        ui->tabWidget->setCurrentIndex(index+1);
    }
    else
    {
        QMessageBox::critical(this, "Try to add more graphs than possible", QString("It is possible to add not more than %0v graphs").arg(NUMBER_OF_PLOTS));
    }
}

void MainWindow::on_PlotButton_clicked()
{
    PlotGraphs();
}


void MainWindow::mousePressEvent(QMouseEvent *event) {
    m_nMouseClick_X_Coordinate = event->x();
    m_nMouseClick_Y_Coordinate = event->y();
}

void MainWindow::mouseMoveEvent(QMouseEvent *event) {
    move(event->globalX()-m_nMouseClick_X_Coordinate,event->globalY()-m_nMouseClick_Y_Coordinate);
}

void MainWindow::on_tabWidget_tabCloseRequested(int index)
{
    QApplication::quit();
}
