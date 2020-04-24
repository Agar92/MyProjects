#ifndef MAINUI_H
#define MAINUI_H
#include <QDialog>

namespace Ui {
  class MainUI;
}

class MainDialog :public QDialog {
  Q_OBJECT
public:
  explicit MainDialog(QWidget *parent = 0);
  ~MainDialog();
public slots:
private:
  Ui::MainUI *m_ui;
};

#endif // MAINUI_H
