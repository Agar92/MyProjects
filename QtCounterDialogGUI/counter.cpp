#include "counter.h"
#include <QPainter>
#include <QTimer>
#include <QDebug>

#include <iostream>
#include "unistd.h"

using namespace std;

Counter::Counter(QWidget *parent)
  :QWidget(parent), m_cap(0), m_timer(0) {
  m_timer = new QTimer(0);
  connect(m_timer, SIGNAL(timeout()), this, SLOT(on_tick()));
}

Counter::~Counter() { delete m_timer; }

void Counter::set(int capacity, int delay) {
  m_timer->stop();

  m_cap = capacity;

  m_val.clear();
  m_val.fill(false, m_cap);

  m_timer->start(delay);
}

void Counter::on_tick() {

/*
//Running red ball - a symbol of process being executed (of loading):
  static int i=m_cap-1;
#ifdef DEBUG
      cout << "1: i:" << i << " i-1:" << i-1 << " i+1:" << i+1 << endl;
      for(int i=0; i<m_cap; ++i) cout << (m_val[i]?1:0) << "  ";
      cout<<endl;
#endif
      m_val[i]=true;
      if(i==m_cap-1) m_val[0]=false;
      else           m_val[i+1]=false;
      if(i==0) i=m_cap; //m_val[1]=false;
      if(i==0) m_val[1]=false;
      --i;
#ifdef DEBUG
      cout << "2: i:" << i << " i-1:" << i-1 << " i+1:" << i+1 << endl;
      for(int i=0; i<m_cap; ++i) cout << (m_val[i]?1:0) << "  ";
      cout<<endl;
      sleep(1);
#endif
*/

//Binary timer:
  //*
  for (int i = m_cap - 1; i >= 0; --i) {
    if (false == m_val[i]) {
      m_val[i] = true;
      break;
    }
    m_val[i] = false;
  }
  //*/
  repaint();
}

void Counter::on_start() { m_timer->start(); }

void Counter::on_stop() { m_timer->stop(); }

/*virtual*/ void Counter::paintEvent(QPaintEvent *) {
  QPainter painter(this);
  if (0 == m_cap) return;

  const int size = qMin(height(), width() / m_cap) - 1;

  for (int i = 0; i < m_cap; ++i) {
    painter.setBrush(QBrush(QColor(255, 0, 0, m_val[i] ? 255 : 0)));
    painter.drawEllipse(i * size, 0, size, size);
  }
}
