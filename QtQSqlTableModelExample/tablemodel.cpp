#include "tablemodel.h"
#include <QColor>

TableModel::TableModel(QObject *parent) : QSqlTableModel(parent)
{

}

QVariant TableModel::data(const QModelIndex &idx, int role) const
{
    if(role == Qt::BackgroundColorRole){
        if(QSqlTableModel::data(this->index(idx.row(), 3)).toInt() < 1.0e3){
            return QColor(Qt::green);
        }
        else if(QSqlTableModel::data(this->index(idx.row(), 3)).toInt() < 1.0e4){
            return QColor(Qt::yellow);
        }
        else if(QSqlTableModel::data(this->index(idx.row(), 3)).toInt() < 1.0e5){
            return QColor(Qt::magenta);
        }
        else if(QSqlTableModel::data(this->index(idx.row(), 3)).toInt() < 1.0e6){
            return QColor(Qt::gray);
        }
        else if(QSqlTableModel::data(this->index(idx.row(), 3)).toInt() < 1.0e7){
            return QColor(Qt::blue);
        }
        else
            return QColor(Qt::red);
    } else if(role == Qt::DisplayRole){
        return QSqlTableModel::data(idx);
    }
    return QVariant();
}
