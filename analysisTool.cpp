#include <iostream>
#include <cmath>
#include <algorithm> 

using namespace std;

int main()
{
    // Constants
    const double pi = 3.14159265358979323846;

    // Input variables
    double ts, f, E, n, fbolt, sOD, Sh, SkH, DeP, CoA_mm, K1, weOfHead, Wins, ti_mm, Weiatt, wP, eWEld, DeOFSk, DeOFSh, b, l;


    cout << "Enter tall vessel thickness (in m): ";
    cin >> ts;

    cout << "Enter f (max allowable stress in MN/m^2): ";
    cin >> f;

    cout << "Enter E (Young's modulus in MN/m^2): ";
    cin >> E;

    cout << "Enter number of bolts: ";
    cin >> n;

    cout << "Enter fbolt (bolt allowable stress in MN/m^2): ";
    cin >> fbolt;

    cout << "Enter shell outside diameter (in m): ";
    cin >> sOD;

    cout << "Enter shell height (in m): ";
    cin >> Sh;

    cout << "Enter skirt height (in m): ";
    cin >> SkH;

    cout << "Enter design pressure (in MN/m^2): ";
    cin >> DeP;

    cout << "Enter corrosion allowance (in mm): ";
    cin >> CoA_mm;

    cout << "Enter K1 (wind load factor): ";
    cin >> K1;

    cout << "Enter weight of head (in kN): ";
    cin >> weOfHead;

    cout << "Enter weight of insulation (in kN): ";
    cin >> Wins;

    cout << "Enter insulation thickness (in mm): ";
    cin >> ti_mm;

    cout << "Enter weight of attachment (in kN): ";
    cin >> Weiatt;

    cout << "Enter wind pressure (in N/m^2): ";
    cin >> wP;

    cout << "Enter weld joint efficiency (0 to 1): ";
    cin >> eWEld;

    cout << "Enter skirt density (in kg/m^3): ";
    cin >> DeOFSk;

    cout << "Enter shell density (in kg/m^3): ";
    cin >> DeOFSh;

    cout << "Enter gusset spacing (in m): ";
    cin >> b;

    cout << "Enter outer radius of bearing plate minus (in m): ";
    cin >> l;

    // Derived values
    double H = Sh + SkH; // Total height

    // Convert corrosion allowance and insulation thickness from mm to meters
    double CoA = CoA_mm / 1000.0;
    double ti = ti_mm / 1000.0;

    // Effective shell diameter including insulation thickness
    double D = sOD + 2 * ti;

    cout << "Effective diameter including insulation D = " << D << " m\n";

    // Calculate inside diameter of shell
    double Di = sOD - 2 * ts;

    // Mean diameter of shell
    double Dm = (sOD + Di) / 2.0;

    // Weight of head in N (from kN)
    double wHead_N = weOfHead * 1000;

    // Weight of shell in N
    double wShell_N = pi * Dm * ts * (H - 4) * 9.81 * DeOFSh;

    cout << "Weight of shell (N): " << wShell_N << endl;

    // Combined weight of shell and head in kN
    double Whead_kN = (wShell_N + wHead_N) / 1000.0;

    cout << "Combined weight of shell and head (kN): " << Whead_kN << endl;

    // Minimum weight (shell + head)
    double Wmin = Whead_kN + wShell_N / 1000.0;

    cout << "Minimum weight Wmin (kN): " << Wmin << endl;

    // Weight of water inside the vessel (assume density 1000 kg/m3)
    double Wwater = pi * (Di / 2) * (Di / 2) * 1000 * 9.81 * Sh / 1000.0;

    cout << "Weight of water (kN): " << Wwater << endl;

    // Maximum weight including attachments and insulation (all in kN)
    double Wmax = Wmin + Weiatt + Wwater + Wins;

    cout << "Maximum weight Wmax (kN): " << Wmax << endl;

    // Choose maximum of Wmin and Wmax for calculations
    double W = max(Wmin, Wmax);

    cout << "Used weight W (kN): " << W << endl;

    // Calculate time factor
    double time = 6.35e-5 * pow(H / D, 1.5) * sqrt(W / ts);

    cout << "Time factor: " << time << endl;

    // Wind load factor K2 = 1 (assumed)
    double K2 = 1.0;

    // Calculate wind force
    double Fbw = K1 * K2 * wP * H * D / 1000.0; // N converted to kN

    cout << "Wind force Fbw (kN): " << Fbw << endl;

    // Wind moment
    double Mw = Fbw * (H / 2);

    cout << "Wind moment Mw (kN.m): " << Mw << endl;

    // Longitudinal bending stress due to wind moment
    double sigWm = 4 * Mw * 7.0 / (pi * Dm * Dm);

    cout << "Wind moment stress sigWm (N/m^2): " << sigWm << endl;

    // Minimum and maximum dead load stresses (kN/m)
    double sigDmin = Wmin / (pi * Dm);

    cout << "Minimum dead load stress sigDmin (kN/m): " << sigDmin << endl;

    double sigDmax = Wmax / (pi * Dm);

    cout << "Maximum dead load stress sigDmax (kN/m): " << sigDmax << endl;

    // Tensile and compressive stresses in MN/m^2 (convert from kN/m^2)
    double sigTen = (sigWm - sigDmin) / 1000.0;
    double sigCom = (sigWm + sigDmax) / 1000.0;

    cout << "Tensile stress sigTen (MN/m^2): " << sigTen << endl;
    cout << "Compressive stress sigCom (MN/m^2): " << sigCom << endl;

    // Calculate z1 and z2 in meters
    double z1 = sigTen / (f * eWEld);
    cout << "z1 (m): " << z1 << endl;

    double z2 = sqrt((sigCom * sOD) / (0.125 * E));
    cout << "z2 (m): " << z2 << endl;

    // Max thickness plus corrosion allowance in meters
    double z = max(z1, z2) + CoA;
    int x = static_cast<int>(z * 1000) + 1; // convert to mm and round up

    cout << "Required shell thickness including corrosion allowance (mm): " << x << endl;

    // Calculate bearing plate area (m^2)
    double area = pi * (sOD - l) * l;
    cout << "Bearing plate area (m^2): " << area << endl;

    // Moment area for gusset
    double za = pi * pow((sOD - l) / 2.0, 2) * l;
    cout << "Moment area za (m^2): " << za << endl;

    // Weight of skirt (in kN)
    double Ws = pi * sOD * (x / 1000.0) * SkH * DeOFSk * 9.81 / 1000.0;
    cout << "Weight of skirt Ws (kN): " << Ws << endl;

    // Bearing compressive stress in MN/m^2
    double sigC = ((Wmax + Ws) / area + (Mw / za)) / 1000.0;
    cout << "Bearing compressive stress sigC (MN/m^2): " << sigC << endl;

    // Thickness of base plate without gussets (m)
    double tbp1 = l * sqrt(3 * sigC / f);
    cout << "Thickness of base plate without gussets (mm): " << tbp1 * 1000 << endl;

    if (tbp1 * 1000 > 20)
    {
        cout << "Gussets are required." << endl;

        // Minimum stress for gusset design
        double sigMin = (Wmin + Ws) / area - (Mw / za);
        cout << "Minimum stress sigMin (kN/m^2): " << sigMin << endl;

        // Moment for gusset
        double My = 0.119 * sigC * l * l;
        cout << "My (moment): " << My << endl;

        // Thickness of base plate with gussets (m)
        double tbp2 = sqrt(6 * My / f);
        cout << "Thickness of base plate with gussets (mm): " << tbp2 * 1000 << endl;

        // Required root area for bolts (m^2)
        double aR = sigMin * area / (n * fbolt * 1000);
        cout << "Required root area aR (m^2): " << aR << endl;

        // Diameter of bolt hole (assuming circular area)
        double dia = sqrt((aR * (-1)) * 28 / 22); // check if negative sign is intentional
        cout << "Bolt hole diameter (mm): " << dia * 1000 << endl;
    }
    else
    {
        cout << "Gussets are not required." << endl;
    }

    return 0;
}

