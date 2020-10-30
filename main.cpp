#include <QCoreApplication>
#include <QtCharts/QLineSeries>
#include <QChart>
#include <QFile>
#include <QTextStream>
#include "eigen-eigen-5a0156e40feb/Eigen/Eigen"

#ifndef M_PI
    #define M_PI 3.14159265358979
#endif
using namespace std;

struct threeDPoint {
    threeDPoint() : x(0), y(0), z(0)
    {

    };
    threeDPoint(double xx, double yy, double zz)
        : x(xx), y(yy), z(zz)
    {
    };
    double x;
    double y;
    double z;
};

void get_plane(const threeDPoint & centroid, const threeDPoint & normal, double & a, double & b, double & c, double & d) {
    d = -(centroid.x * normal.x + centroid.y * normal.y + centroid.z * normal.z);
    a = normal.x;
    b = normal.y;
    c = normal.z;
    //std::cout << "The plane is " << weighted_dir.x << "*x + " << weighted_dir.y << "*y + " << weighted_dir.z << "*z + " << d << "=0" << std::endl;
}

double dot_product(threeDPoint a, threeDPoint b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void normalize(threeDPoint & vec) {
    double len = vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
    len = sqrt(len);
    if (len <= 0)
        return;
    vec.x /= len;
    vec.y /= len;
    vec.z /= len;
}

double cal_error(const std::vector<threeDPoint> & points, double a, double b, double c, double d) {
    double error = 0;
    double denum = a * a + b * b + c * c;
    if (denum <  0.0001)
        return error;
    for (auto it = points.begin(); it != points.end(); ++it) {
        double numerator = std::fabs(a * it->x + b * it->y + c * it->z + d);
        error += numerator / denum;
    }
    return error;
}

threeDPoint planeFitting(std::vector<threeDPoint> points) {
    double xx_sum = 0, yy_sum = 0, zz_sum = 0;
    int total = (int)points.size();
    for (int i = 0; i < total; ++i) {
        xx_sum += points[i].x;
        yy_sum += points[i].y;
        zz_sum += points[i].z;
    }
    double cx = xx_sum / total;
    double cy = yy_sum / total;
    double cz = zz_sum / total;
    threeDPoint centroid(cx, cy, cz);
    double xx = 0, xy = 0, xz = 0;
    double yy = 0, yz = 0, zz = 0;
    for (int i = 0; i < total; ++i) {
        double dev_x = points[i].x - cx;
        double dev_y = points[i].y - cy;
        double dev_z = points[i].z - cz;
        xx += dev_x * dev_x;
        xy += dev_x * dev_y;
        xz += dev_x * dev_z;
        yy += dev_y * dev_y;
        yz += dev_y * dev_z;
        zz += dev_z * dev_z;
    }
    xx /= total;
    xy /= total;
    xz /= total;
    yy /= total;
    yz /= total;
    zz /= total;
    // Method 1
    threeDPoint weighted_dir = threeDPoint(0, 0, 0);
    {
        double det_x = yy * zz - yz * yz;
        double axis_x = det_x;
        double axis_y = xz * yz - xy * zz;
        double axis_z = xy * yz - xz * yy;
        threeDPoint axis_dir(axis_x, axis_y, axis_z);
        //std::cout << "axis_dir: (" << axis_dir.x << ", " << axis_dir.y << ", " << axis_dir.z << ")" << std::endl;
        double weight = det_x * det_x;
        if (dot_product(weighted_dir, axis_dir) < 0.0) {
            weight = -weight;
        }
        weighted_dir.x += axis_dir.x * weight;
        weighted_dir.y += axis_dir.y * weight;
        weighted_dir.z += axis_dir.z * weight;
    }
    {
        double det_y = xx * zz - xz * xz;
        double axis_x = xz * yz - xy * zz;
        double axis_y = det_y;
        double axis_z = xy * xz - yz * xx;
        threeDPoint axis_dir(axis_x, axis_y, axis_z);
        //std::cout << "axis_dir: (" << axis_dir.x << ", " << axis_dir.y << ", " << axis_dir.z << ")" << std::endl;
        double weight = det_y * det_y;
        if (dot_product(weighted_dir, axis_dir) < 0.0) {
            weight = -weight;
        }
        weighted_dir.x += axis_dir.x * weight;
        weighted_dir.y += axis_dir.y * weight;
        weighted_dir.z += axis_dir.z * weight;
    }
    {
        double det_z = xx * yy - xy * xy;
        double axis_x = xy * yz - xz * yy;
        double axis_y = xy * xz - yz * xx;
        double axis_z = det_z;
        threeDPoint axis_dir(axis_x, axis_y, axis_z);
        //std::cout << "axis_dir: (" << axis_dir.x << ", " << axis_dir.y << ", " << axis_dir.z << ")" << std::endl;
        double weight = det_z * det_z;
        if (dot_product(weighted_dir, axis_dir) < 0.0) {
            weight = -weight;
        }
        weighted_dir.x += axis_dir.x * weight;
        weighted_dir.y += axis_dir.y * weight;
        weighted_dir.z += axis_dir.z * weight;
        normalize(weighted_dir);
        double a, b, c, d;
        get_plane(centroid, weighted_dir, a, b, c, d);
        //double error = cal_error(points, a, b, c, d);
        //Convert to degree
        weighted_dir.x *= 180.0/M_PI;
        weighted_dir.y *= 180.0/M_PI;

        qInfo("The plane is %f*x + %f*y + %f*z + %f*d = 0", a, b, c, d);
        {
            qInfo("The plane is %f*x + %f*y + %f*d = z", a/-c, b/-c, d/-c);
        }
        //DBOUT("The plane is " << a << "*x + " << b << "*y + " << c << "*z + " << d << "=0" << std::endl);
        //DBOUT("Normal vector: (" << weighted_dir.x << ", " << weighted_dir.y << ", " << weighted_dir.z << ")" << std::endl);
        //qInfo("The plane is weighted_dir.x: %f weighted_dir.y: %f weighted_dir.z: %f", weighted_dir.x, weighted_dir.y, weighted_dir.z);
        //qInfo("The plane error: %f", error);
        //double error = cal_error(points, a, b, c, d);
        //DBOUT( "Total error: " << error << std::endl );
        return weighted_dir;
    }
}

vector<double> fitCurve(const vector<double> & x, const vector<double> & y, int order)
{
    size_t n = x.size();
    double minX = 999999;
    double maxX = -999999;
    for (size_t i = 0; i < n; i++) {
        minX = std::min(minX, x[i]);
        maxX = std::max(maxX, x[i]);
    }
    Eigen::MatrixXd X(n, order + 1);
    for (size_t i = 0; i < n; ++i) {
        double tmp = 1;
        //order次方
        for (int j = 0; j <= order; ++j) {
            X(i, j) = tmp;
            tmp *= x[i];
        }
    }
    Eigen::MatrixXd Y(n, 1);
    for (size_t i = 0; i < n; ++i) {
        Y(i, 0) = y[i];
    }
    Eigen::MatrixXd Xt(order + 1, n);
    Xt = X.transpose();
    Eigen::MatrixXd XtX(order + 1, order + 1);
    XtX = Xt * X;
    Eigen::MatrixXd XtX_inverse(order + 1, order + 1);
    XtX_inverse = XtX.inverse();
    Eigen::MatrixXd A(order + 1, 1);
    A = XtX_inverse * Xt * Y;
    Eigen::MatrixXd B(order + 1, 1);
    B = X * A;
    Eigen::MatrixXd Ans(n, 1);
    Ans = X * A;
    double error = 0;
    double average = 0;
    for (size_t i = 0; i < n; ++i) {
        average += (B(i, 0) - Y(i, 0));
    }
    average = average/n;
    for (size_t i = 0; i < n; ++i) {
        double diff = B(i, 0) - Y(i, 0) - average;
        error +=(diff * diff);
    }
    vector<double> ans;
    for (int i = 0; i <= order; ++i) {
        ans.push_back(A(i, 0));
    }
    return ans;
}

void curveFittingUnitTest()
{
    double freq_scale = 500;

    QtCharts::QLineSeries series;

    QFile file("C:\\Users\\emil\\Desktop\\1.csv");

    if (file.open(QIODevice::ReadOnly))
    {
        QTextStream in(&file);
        int length = 0;
        while(!in.atEnd()){
            QString line = in.readLine();
            length++;
            if (length > 1) { //Skip header line
                QStringList list = line.split(',');
                if (list.size() == 2) {
                    double x = list[0].toDouble();
                    double y = list[1].toDouble();
                    series.append(x, y);
                }
            }
        }
        file.close();
    }

    vector<double> x, y;
    QVector<QPointF> pts(series.pointsVector());
    for (int pi=0; pi < pts.size(); pi++) {
        pts[pi] = QPointF(pts[pi].x(), pts[pi].y());
        qInfo("x: %f y: %f", pts[pi].x(), pts[pi].y());
        x.push_back( pts[pi].x() * freq_scale);
        y.push_back( pts[pi].y());
    }
    vector<double> coefficient = fitCurve(x, y, 5);

    //Unit test
    double x_point = 125;
    double y_point = 0;
    for (int i = 0; i < coefficient.size(); i++) {
        qInfo("c: %d %f", i, coefficient[i]);
        y_point += pow(x_point, i)*coefficient[i];
    }
    qInfo("x: %f y: %f", x_point, y_point);
}

void planeFittingUnitTest()
{
    double x[10] = { -1.3665096083681405, 0.932377583174361, -3.484999901434442, -3.734641761539269, 3.62312454263758, 3.5958770103915754, -1.7253324302895736, -1.319650743721322, -3.0910801380675537, 1.2055393876976073};
    double y[10] = { -1.753677437590138, 0.2357532087900518, 1.476233648961962, -0.3802579066654683, -2.2454264290377166, -2.187759869190746, -2.1088223752647197, 3.5361076037682846, -1.796371681391399, 3.717959144424345};
    double z[10] = { -2.994051529506695, 7.5720147927188775, 2.4587011440170023, -3.6100572430749427, 5.50996979816201, 5.628474413210913, -4.7771319863733055, 12.96902132386221, -6.571275320309304, 18.56495620866825 };

    double x1[10] = {4.741494561013635, 1.3455734190842508, -2.3873836114654843, 2.237083602376602, -2.0852357734780353, -2.479139647087173, -2.916868948078637, -1.8282375104853248, -3.4394712225601225, 1.263280302841233};
    double y1[10] = {3.5706442738418147, -3.0656122362760794, -2.568068810740603, -2.098898473235116, 0.5628788713827468, 4.861633577635269, -2.281968082741841, -3.563364703711649, 2.9047541079064434, 3.396824422488935};
    double z1[10] = {23.922433729452763, 1.344350868044597, -5.882672624398617, 2.3181296169220675, 3.9517457623203383, 15.274051909657869, -8.566948363065464, -9.010840739989446, 6.880479921801626, 18.223370989944343};

    qInfo("Test Set 1");
    {
        vector<threeDPoint> points;
        for (int i = 0; i < 10; i++){
            points.push_back(threeDPoint(x[i], y[i], z[i]));
        }
        threeDPoint weighted_vector = planeFitting(points);
        double x_tilt = weighted_vector.z*weighted_vector.x;
        double y_tilt = weighted_vector.z*weighted_vector.y;
        qInfo("xTilt: %f yTilt: %f", x_tilt, y_tilt);
    }

    qInfo("Test Set 2");
    {
        vector<threeDPoint> points;
        for (int i = 0; i < 10; i++){
            points.push_back(threeDPoint(x1[i], y1[i], z1[i]));
        }
        threeDPoint weighted_vector = planeFitting(points);
        double x_tilt = weighted_vector.z*weighted_vector.x;
        double y_tilt = weighted_vector.z*weighted_vector.y;
        qInfo("xTilt: %f yTilt: %f", x_tilt, y_tilt);
    }
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    //Curve Fitting
    //curveFittingUnitTest();

    //Plane Fitting
    planeFittingUnitTest();

    return a.exec();
}
