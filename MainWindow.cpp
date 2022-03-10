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
    QLayout* layout = new QHBoxLayout;
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
    connect(ui->pgsRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setSolver(kPGS); });
    connect(ui->bppRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setSolver(kBPP); });

    connect(ui->boxPlane, SIGNAL(clicked()), viewer, SLOT(createBoxOnPlane()));
    connect(ui->boxInclinedPlane, SIGNAL(clicked()), viewer, SLOT(createBoxOnInclinedPlane()));
    connect(ui->boxBallStack, SIGNAL(clicked()), viewer, SLOT(createBoxBallStack()));
    connect(ui->bunnies, SIGNAL(clicked()), viewer, SLOT(createBunnies()));
    connect(ui->marblesBox, SIGNAL(clicked()), viewer, SLOT(createMarbleBox()));

    connect(ui->cornersRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setNormalSamplingMethod(kCorners); });
    connect(ui->gridRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setNormalSamplingMethod(kGrid); });
    connect(ui->randomRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setNormalSamplingMethod(kRandom); });

    connect(ui->noneRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionSamplingMethod(kNone); });
    connect(ui->aggCoeffRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionSamplingMethod(kAggCoeff); });
    connect(ui->ForceEquivRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionSamplingMethod(kForceEquiv); });

    connect(ui->uniformRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionDistribution(kUniform); });
    connect(ui->gaussianRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionDistribution(kGaussian); });
    connect(ui->rampRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionDistribution(kRamp); });
    connect(ui->stepRadio, &QRadioButton::toggled, viewer, [viewer] { viewer->setFrictionDistribution(kStep); });


    // Update status bar message
    connect(viewer, SIGNAL(statusMessageChanged(QString)), ui->statusBar, SLOT(showMessage(QString)));
}