#include "settingsform.h"
#include "ui_settingsform.h"

#include <iostream>
#include <sstream>

SettingsForm::SettingsForm(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SettingsForm)
{
    ui->setupUi(this);
    connect(ui->ApplySettingsButton, SIGNAL(clicked()), this, SLOT(on_ApplySettingsButton_clicked()));
}

SettingsForm::~SettingsForm()
{
    delete ui;
}

QVector<double> SettingsForm::GetSettings()
{
    std::cout<<"SettingsForm::GetSettings() BEGIN:"<<std::endl;
    settings.clear();
    std::stringstream s1(ui->lineEditXmin->text().toUtf8().constData());
    double xmin;
    s1>>xmin;
    std::stringstream s2(ui->lineEditXmax->text().toUtf8().constData());
    double xmax;
    s2>>xmax;
    std::stringstream s3(ui->lineEditYmin->text().toUtf8().constData());
    double ymin;
    s3>>ymin;
    std::stringstream s4(ui->lineEditYmax->text().toUtf8().constData());
    double ymax;
    s4>>ymax;
    settings.push_back(xmin);
    settings.push_back(xmax);
    settings.push_back(ymin);
    settings.push_back(ymax);
    return settings;
}

void SettingsForm::on_ApplySettingsButton_clicked()
{
    QVector<double> settings = GetSettings();
    std::cout<<"emit SettingsChanged():"<<std::endl;
    emit SettingsChanged(settings);
}
