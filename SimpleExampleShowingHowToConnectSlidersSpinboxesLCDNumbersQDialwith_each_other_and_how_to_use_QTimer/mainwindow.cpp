#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QSqlQueryModel>
#include <QtSql/QSqlQuery>

#include <QDebug>
#include <QMessageBox>
#include <QSqlRecord>
#include <QtSql/QSqlError>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

    qDebug() << "STEP #1";

    ui->setupUi(this);
    QSqlQueryModel * model = new QSqlQueryModel();
    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE");
        db.setDatabaseName("mydb.sqlite");
        if (!db.open()) {
            qDebug() << "Database not opened!";
        }

    qDebug() << "STEP #2";

    QSqlQuery query(db);
    query.exec("select * from one;");
    QSqlRecord rec = query.record();
    int id, age;
    QString name;
    qDebug() << query.next();
    while (query.next()) {
        id = query.value(rec.indexOf("id")).toInt();
        name = query.value(rec.indexOf("name")).toString();
        age = query.value(rec.indexOf("address")).toInt();

        qDebug() << "id is " << id
                 << ". name is " << name
                 << ". age is" << age
                 << query.lastError().text()
                 << query.lastError().number();
        }

    qDebug() << "STEP #3";

    model->setQuery(query);
    ui->tableView->setModel(model);

    qDebug() << "STEP #4";


    ui->progressBar->setMinimum(0);
    ui->progressBar->setMaximum(100);
    ui->progressBar->setValue(0);
    ui->horizontalSlider->setMinimum(0);
    ui->horizontalSlider->setMaximum(100);
    ui->horizontalSlider->setValue(0);
    ui->spinBox->setMinimum(0);
    ui->spinBox->setMaximum(100);
    ui->spinBox->setValue(0);
    ui->lcdNumber->display(0);

    //ConnectSignals();

    timer = new QTimer();
    timer->setInterval(10);//Задаем интервал таймера
    //connect(timer, SIGNAL(timeout()), this, SLOT(UpdateTime())); // Подключаем сигнал таймера к нашему слоту
    timer->start();

    connect(timer, SIGNAL(timeout()), this, SLOT(UpdateTime()));

    ui->statusbar->showMessage("Hello!", 3000);

    QTimer::singleShot(5000, ui->lcdNumber, SLOT(hide()));

    QTimer::singleShot(10000, ui->lcdNumber, SLOT(show()));

    ui->dial->setRange(0, 360);
    ui->dial->setNotchTarget(20);
    ui->dial->setNotchesVisible(true);
    ui->dial->setWrapping(true);


}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::ConnectSignals()
{
    connect(ui->progressBar, &QProgressBar::valueChanged, this, &MainWindow::Valuechanged);
    connect(ui->horizontalSlider, &QSlider::valueChanged, this, &MainWindow::Valuechanged);
    connect(ui->spinBox, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, &MainWindow::Valuechanged);
    connect(ui->doubleSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &MainWindow::Valuechanged);
    connect(ui->verticalSlider, &QSlider::valueChanged, this, &MainWindow::Valuechanged);
}

void MainWindow::DisconnectSignals()
{
    disconnect(ui->progressBar, &QProgressBar::valueChanged, this, &MainWindow::Valuechanged);
    disconnect(ui->horizontalSlider, &QSlider::valueChanged, this, &MainWindow::Valuechanged);
    disconnect(ui->spinBox, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, &MainWindow::Valuechanged);
    disconnect(ui->doubleSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, &MainWindow::Valuechanged);
    disconnect(ui->verticalSlider, &QSlider::valueChanged, this, &MainWindow::Valuechanged);
}

void MainWindow::on_ViewButton_clicked()
{

    QMessageBox::information(this, "View", "View button clicked");
}

void MainWindow::Valuechanged(int val)
{
    static int i=0;
    //DisconnectSignals();
    ui->progressBar->setValue(val);
    ui->horizontalSlider->setValue(val);
    ui->spinBox->setValue(val);
    ui->doubleSpinBox->setValue((double)val);
    ui->lcdNumber->display(val);
    ui->verticalSlider->setValue(val);
    //ConnectSignals();
    ++i;
}

void MainWindow::UpdateTime()
{
    static int i=0;
    static bool b=true;
    qDebug() << "UpdateTime():" ;

    //ui->horizontalSlider->setValue(i);
    Valuechanged(i);


    if(i==99) b=false;
    if(i==0)   b=true;
    if(b) ++i;
    else  --i;
}

void MainWindow::on_dial_valueChanged(int value)
{
    ui->lcdNumber_2->display(value);
}
