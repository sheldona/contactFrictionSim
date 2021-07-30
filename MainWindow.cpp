#include "MainWindow.h"
#include "ui_mainwindow.h"

#include "viewer/SimViewer.h"
#include <QBoxLayout>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::addViewer(SimViewer* viewer)
{
    // Set the viewer
    QLayout *layout = new QHBoxLayout;
    layout->addWidget(viewer);
    ui->frame->setLayout(layout);

    connect(ui->stepOnce, SIGNAL(clicked()), viewer, SLOT(stepOnce()));
    connect(ui->maxIterations, SIGNAL(valueChanged(int)), viewer, SLOT(setMaxIterations(int)));
    connect(ui->drawMesh, SIGNAL(toggled(bool)), viewer, SLOT(setDrawMesh(bool)));
    connect(ui->drawContacts, SIGNAL(toggled(bool)), viewer, SLOT(setDrawContacts(bool)));
    connect(ui->timeStep, SIGNAL(valueChanged(double)), viewer, SLOT(setTimestep(double)));
    connect(ui->subSteps, SIGNAL(valueChanged(int)), viewer, SLOT(setSubsteps(int)));
    connect(ui->contactStiffness, SIGNAL(valueChanged(double)), viewer, SLOT(setContactStiffness(double)));
    connect(ui->contactDamping, SIGNAL(valueChanged(double)), viewer, SLOT(setContactDamping(double)));
    connect(ui->frictionCoeff, SIGNAL(valueChanged(double)), viewer, SLOT(setFrictionCoefficient(double)));
    connect(ui->pause, SIGNAL(toggled(bool)), viewer, SLOT(setPaused(bool)));

    connect(ui->boxPlane, SIGNAL(clicked()), viewer, SLOT(createBoxOnPlane()));
    connect(ui->boxBallStack, SIGNAL(clicked()), viewer, SLOT(createBoxBallStack()));
    connect(ui->bunnies, SIGNAL(clicked()), viewer, SLOT(createBunnies()));
    connect(ui->marblesBox, SIGNAL(clicked()), viewer, SLOT(createMarbleBox()));

    // Update status bar message
    connect(viewer, SIGNAL(statusMessageChanged(QString)), ui->statusBar, SLOT(showMessage(QString)));
}
