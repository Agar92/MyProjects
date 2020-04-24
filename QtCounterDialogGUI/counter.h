#ifndef COUNTER_H
#define COUNTER_H
#include <QWidget>
#include <QVector>

//!< двоичный счетчик
class Counter :public QWidget {
  Q_OBJECT
public:
  Counter(QWidget* parent);
  virtual ~Counter();
  void set(int capacity, int delay);
              //!< установка разрядности и задержки таймера
protected slots:
  void on_tick();
              //!< обработка сигнала таймера
  void on_start();
              //!< обработка сигнала запуска
  void on_stop();
              //!< обработка сигнала остановки
protected:
  virtual void paintEvent(QPaintEvent *);

  int m_cap;
              //!< количество разрядов
  QVector<bool> m_val;
              //!< отображаемое значение
  QTimer *m_timer;
              //!< таймер
};

#endif // COUNTER_H
