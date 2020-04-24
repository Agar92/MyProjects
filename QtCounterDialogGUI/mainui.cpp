#include "mainui.h"
#include "ui_mainui.h"

MainDialog::MainDialog(QWidget *parent)
  : QDialog(parent), m_ui(new Ui::MainUI) { 
  m_ui->setupUi(this);

  m_ui->m_counter->set(8, 200);

  connect(m_ui->m_start, SIGNAL(clicked()), m_ui->m_counter, SLOT(on_start()));
  connect(m_ui->m_stop, SIGNAL(clicked()), m_ui->m_counter, SLOT(on_stop()));
}

MainDialog::~MainDialog() {
  delete m_ui;
}
