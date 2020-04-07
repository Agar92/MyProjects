#ifndef SETTINGSFORM_H
#define SETTINGSFORM_H

#include <QWidget>

namespace Ui {
class SettingsForm;
}

class SettingsForm : public QWidget
{
    Q_OBJECT

public:
    explicit SettingsForm(QWidget *parent = nullptr);
    ~SettingsForm();
    QVector<double> GetSettings();



signals:
    void SettingsChanged(QVector<double> &);

private slots:
    void on_ApplySettingsButton_clicked();

private:
    Ui::SettingsForm *ui;
    QVector<double> settings;
};

#endif // SETTINGSFORM_H
