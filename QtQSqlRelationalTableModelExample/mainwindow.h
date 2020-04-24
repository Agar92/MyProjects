#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSqlRelationalTableModel>
#include <QSqlRelationalDelegate>
#include <QSqlRelation>

/* Подлючаем заголовчный файл для работы с базой данных */
#include "database.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow  *ui;
    /* В проекте используются объекты для работы с базой данных
     * и моделью представления таблицы базы данных
     * */
    DataBase                    *db;
    QSqlRelationalTableModel    *modelMain;
    QSqlRelationalTableModel    *modelDevice;

private:
    /* Также присутствуют два метода, которые формируют модель
     * и внешний вид TableView
     * */
    void setupMainModel(const QString &tableName, const QStringList &headers);
    void setupDeviceModel(const QString &tableName, const QStringList &headers);
    void createUI();
};

#endif // MAINWINDOW_H
