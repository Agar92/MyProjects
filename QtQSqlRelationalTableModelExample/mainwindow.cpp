#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    /* Первым делом необходимо создать объект для работы с базой данных
     * и инициализировать подключение к базе данных
     * */
    db = new DataBase();
    db->connectToDataBase();

    /* После чего производим наполнение таблицы базы данных
     * контентом, который будет отображаться в tableView и tableViewDevice
     * */
    for(int i = 1; i < 4; i++){
        QVariantList data;
        data.append("Device " + QString::number(i));
        data.append("192.168.0." + QString::number(i));
        db->insertIntoDeviceTable(data);
    }

    for(int i = 0; i < 10; i++){
        QVariantList data;
        QString random = QString::number(qrand() % ((4 + 1) - 1) + 1);
        data.append(QDate::currentDate());
        data.append(QTime::currentTime());
        data.append(random);
        data.append(random);
        db->insertIntoMainTable(data);
    }

    /* Инициализируем модели для представления данных
     * с заданием названий колонок
     * */
    this->setupMainModel(TABLE,
                     QStringList() << trUtf8("id")
                                   << trUtf8("Дата")
                                   << trUtf8("Время")
                                   << trUtf8("Имя хоста")
                                   << trUtf8("IP адрес")
               );

    this->setupDeviceModel(DEVICE,
                     QStringList() << trUtf8("id")
                                   << trUtf8("Имя хоста")
                                   << trUtf8("IP адрес")
               );
    /* Инициализируем внешний вид таблицы с данными
     * */
    this->createUI();
}

MainWindow::~MainWindow()
{
    delete ui;
}

/* Метод для инициализации модели представления данных
 * */
void MainWindow::setupMainModel(const QString &tableName, const QStringList &headers)
{
    /* Производим инициализацию модели представления данных
     * с установкой имени таблицы в базе данных, по которому
     * будет производится обращение в таблице
     * */
    modelMain = new QSqlRelationalTableModel(this);
    modelMain->setTable(tableName);
    /* Устанавливаем связи с таблицей устройств, по которым будет производится
     * подстановка данных
     * В метода setRelation указывается номер колонки, в которой будет
     * производится подстановка, а также с помощью класса
     * QSqlRelation указывается имя таблицы,
     * параметр, по которому будет произведена выборка строки
     * и колонка, из которой будут взяты данные
     * */
    modelMain->setRelation(3, QSqlRelation(DEVICE, "id", DEVICE_HOSTNAME));
    modelMain->setRelation(4, QSqlRelation(DEVICE, "id", DEVICE_IP));

    /* Устанавливаем названия колонок в таблице с сортировкой данных
     * */
    for(int i = 0; i < modelMain->columnCount(); i++){
        modelMain->setHeaderData(i,Qt::Horizontal,headers[i]);
    }
    // Устанавливаем сортировку по возрастанию данных по нулевой колонке
    modelMain->setSort(0,Qt::AscendingOrder);
    modelMain->select(); // Делаем выборку данных из таблицы
}

void MainWindow::setupDeviceModel(const QString &tableName, const QStringList &headers)
{
    /* Производим инициализацию модели представления данных
     * с установкой имени таблицы в базе данных, по которому
     * будет производится обращение в таблице
     * */
    modelDevice = new QSqlRelationalTableModel(this);
    modelDevice->setTable(tableName);

    /* Устанавливаем названия колонок в таблице с сортировкой данных
     * */
    for(int i = 0; i < modelDevice->columnCount(); i++){
        modelDevice->setHeaderData(i,Qt::Horizontal,headers[i]);
    }
    // Устанавливаем сортировку по возрастанию данных по нулевой колонке
    modelDevice->setSort(0,Qt::AscendingOrder);
    modelDevice->select(); // Делаем выборку данных из таблицы
}

void MainWindow::createUI()
{
    ui->tableView->setModel(modelMain);     // Устанавливаем модель на TableView
    ui->tableView->setColumnHidden(0, true);    // Скрываем колонку с id записей
    // Разрешаем выделение строк
    ui->tableView->setSelectionBehavior(QAbstractItemView::SelectRows);
    // Устанавливаем режим выделения лишь одно строки в таблице
    ui->tableView->setSelectionMode(QAbstractItemView::SingleSelection);
    // Устанавливаем размер колонок по содержимому
    ui->tableView->resizeColumnsToContents();
    ui->tableView->setItemDelegate(new QSqlRelationalDelegate(ui->tableView));
    //ui->tableView->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->tableView->horizontalHeader()->setStretchLastSection(true);

    modelMain->select(); // Делаем выборку данных из таблицы

    ui->tableViewDevice->setModel(modelDevice);     // Устанавливаем модель на TableView
    ui->tableViewDevice->setColumnHidden(0, true);    // Скрываем колонку с id записей
    // Разрешаем выделение строк
    ui->tableViewDevice->setSelectionBehavior(QAbstractItemView::SelectRows);
    // Устанавливаем режим выделения лишь одно строки в таблице
    ui->tableViewDevice->setSelectionMode(QAbstractItemView::SingleSelection);
    // Устанавливаем размер колонок по содержимому
    ui->tableViewDevice->resizeColumnsToContents();
    ui->tableViewDevice->setItemDelegate(new QSqlRelationalDelegate(ui->tableViewDevice));
    //ui->tableViewDevice->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->tableViewDevice->horizontalHeader()->setStretchLastSection(true);

    modelDevice->select();
}
