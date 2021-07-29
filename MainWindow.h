#pragma once

/**
 * @file MainWindow.h
 *
 * @brief Main Qt Window.
 *
 */

#include <QMainWindow>

QT_FORWARD_DECLARE_CLASS(SimViewer)

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void addViewer(SimViewer* viewer);

private:
    Ui::MainWindow *ui;
};
