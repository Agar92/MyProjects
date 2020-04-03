#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "parser.h"
#include "qcustomplot.h"

#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->centralwidget->setMinimumSize(QSize(400, 600));

    ui->progressBar->setValue(0);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_PlotButton_clicked()
{

    tokens texpr, pexpr;
    Variables expvars;
    Massives expmasvars;
    string expr = ui->lineEdit->text().toUtf8().constData();;

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

      const int N=10000;
      QVector<double> x(N), y(N);
      for(int i=0; i<N; ++i) x[i]=0.1+i * 100.0/(N-1);
      for(int i=0; i<N; ++i)
      {
          expvars["x"]=x[i];
          cout << "Result:" << endl;
          y[i]=ResultExpr(pexpr, expvars, expmasvars);
          cout << y.at(i) << endl;
          if(i%100==0)
          {
              ui->progressBar->setValue((int)i/100 + 1);
          }
      }

      ui->widget->clearGraphs();//Если нужно, но очищаем все графики
      //Добавляем один график в widget
      ui->widget->addGraph();
      //Говорим, что отрисовать нужно график по нашим двум массивам x и y
      ui->widget->graph(0)->setData(x, y);

      if(ui->lineEdit2->text().size()>0)
      {
        tokens texpr2, pexpr2;
        Variables expvars2;
        Massives expmasvars2;
        string expr2 = ui->lineEdit2->text().toUtf8().constData();;
        CreateTokensFromExpression(expr2, texpr2);
        CreatePostfixFromTokens(texpr2, pexpr2);
        QVector<double> y1(N, 1);
        //for(int i=0; i<N; ++i) y1[i]=x[i];
        for(int i=0; i<N; ++i)
        {
            expvars2["x"]=x[i];
            cout << "Result:" << endl;
            y1[i]=ResultExpr(pexpr2, expvars2, expmasvars2);
            cout << y1.at(i) << endl;
        }

        //Добавляем один график в widget
        ui->widget->addGraph();
        //Говорим, что отрисовать нужно график по нашим двум массивам x и y
        ui->widget->graph(1)->setData(x, y1);
        QPen pen;  // создаёт перо по умолчанию
        pen.setStyle(Qt::DashDotLine);
        pen.setWidth(3);
        pen.setBrush(Qt::green);
        ui->widget->graph(1)->setPen(pen);
      }

      //Подписываем оси Ox и Oy
      ui->widget->xAxis->setLabel("x");
      ui->widget->yAxis->setLabel("y");
      double minX = x[0], maxX = x[N-1];
      //Установим область, которая будет показываться на графике
      ui->widget->xAxis->setRange(minX, maxX);//Для оси Ox
      //Для показа границ по оси Oy сложнее, так как надо по правильному
      //вычислить минимальное и максимальное значение в векторах
      double minY = y[0], maxY = y[0];
      for (int i=1; i<N; i++)
      {
          if (y[i]<minY) minY = y[i];
          if (y[i]>maxY) maxY = y[i];
      }
      ui->widget->yAxis->setRange(minY, maxY);//Для оси Oy
      //И перерисуем график на нашем widget
      ui->widget->replot();
      //QTimer::singleShot(2000, this, SLOT(ui->progressBar->setValue(0)));
      int I=0;
      for(int i=0; i<1000000000; ++i) I+=i;
      ui->progressBar->setValue(0);
    }
  }
