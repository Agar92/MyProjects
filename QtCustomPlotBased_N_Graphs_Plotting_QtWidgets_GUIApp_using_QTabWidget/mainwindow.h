#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QString>
#include <QVector>
#include <QMap>
#include "addgraphform.h"
#include "settingsform.h"

#include "parser.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    Ui::MainWindow * GetMainWindow(){return ui;}

    void PlotGraphs();

private slots:
    void on_PlotButton_clicked();
    void on_AddGraphButton_clicked();
    void on_tabWidget_Close();

    void on_tabWidget_tabCloseRequested(int index);

public slots:
    void ApplySettingChanges(QVector<double> & vec);


private:
    static const int NUMBER_OF_PLOTS=3;
    int AddGraphTabsNumber;
    Ui::MainWindow *ui;
    AddGraphForm *form[NUMBER_OF_PLOTS];
    SettingsForm * settingsForm;
    bool alive[NUMBER_OF_PLOTS];
    QColor colors[NUMBER_OF_PLOTS];
    QColor COLORS_STORED[9]{QColor(Qt::black), QColor(Qt::red), QColor(Qt::cyan),
                            QColor(Qt::green), QColor(Qt::blue), QColor(Qt::magenta),
                            QColor(Qt::yellow), QColor(Qt::darkGray), QColor(Qt::darkRed)};
    int ColorCounter;
    int CreatedGraphsCounter;
    //plot ranges:
    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;

    //this is for moving the window with
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    int m_nMouseClick_X_Coordinate;
    int m_nMouseClick_Y_Coordinate;

};
#endif // MAINWINDOW_H
