#include <QApplication>

#include "MainWindow.h"
#include "viewer/SimViewer.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    // Set the core profile and version of OpenGL shaders.
    //
    QSurfaceFormat fmt;
    fmt.setVersion(4, 0);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    fmt.setDepthBufferSize(32);
    fmt.setSamples(0);
    QSurfaceFormat::setDefaultFormat(fmt);

    MainWindow w;

    // Instantiate the simulation viewer.
    //
    SimViewer v;
    w.addViewer(&v);
    w.show();

    return a.exec();
}
