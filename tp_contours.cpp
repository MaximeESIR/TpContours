/****************************************************************************
 * Copyright (C) 2016 Universite de Rennes 1. All rights reserved.
 *
 * This software was developed at:
 * Universite de Rennes 1
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 *
 * This file uses the ViSP library.
 *****************************************************************************/

/****************************************************************************
 * NOMS - PRENOMS:
 *  -Alric Théo
 *	-Blaret Maxime
 *
 * Date :
 *****************************************************************************/

#include <iostream>

#include <visp/vpConfig.h>
#include <visp/vpDebug.h>

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpDisplayX.h>

using namespace std;

/**
 * @brief affiche une image � l'�cran � la position (posX,posY) et attend un clic
 * @param img : l'image � afficher
 * @param posX, posY : coordonn�es spatiales pour positionner affichage de la fenetre sur l'�cran
 * @param title : titre de la fenetre graphique
 * @warning : fonction bloquante
 */
void afficheImage(vpImage<unsigned char> img, int posX, int posY, const char *title)
{
    vpDisplayX d(img, posX, posY, title);
    vpDisplay::display(img);
    vpDisplay::flush(img);
    vpDisplay::getClick(img);
    vpDisplay::close(img);
}
// Fonction d'affichage d'une image sur sortie standard
void affiche(const vpImage<unsigned char> &src)
{
    int i, j;
    printf("\n ");
    for (i = 0; i < src.getHeight(); i++)
    {
        for (j = 0; j < src.getWidth(); j++)
        {
            printf("%d ", src[i][j]);
        }
        printf("\n ");
    }
}

void filtrage2D(const vpImage<unsigned char> &I, vpImage<double> &Ic, const vpMatrix &K)
{
    // Allocate memory for an [height x width] image and initialize the image to val.
    Ic.resize(I.getHeight(), I.getWidth(), 0);
    int offset = K.getRows() / 2; // taille du filtre /2
    double val = 0;
    for (int i = 0; i < I.getHeight(); i++)
    { // Parcours de l'image
        for (int j = 0; j < I.getWidth(); j++)
        {
            val = 0; // Valeur de I'[i,j]
            for (int s = 0; s < K.getRows(); s++)
            { // Parcours du filtre
                for (int t = 0; t < K.getRows(); t++)
                {
                    if (i - offset + s < 0 || j - offset + t < 0 || i - offset + s >= I.getHeight() || j - offset + t >= I.getWidth()) // Pour éviter pb de segmentation
                        val += 0;
                    else
                        val += (K[s][t]) * I[i - offset + s][j - offset + t];
                }
            }
            if (val > 255)
            { // Eviter de corrompre l'image
                Ic[i][j] = 255;
            }
            else
            {
                Ic[i][j] = val;
            }
        }
    }
}
// Calcul du gradient
void gradient(vpImage<unsigned char> I0, vpImage<double> &IGradX, vpImage<double> &IGradY, vpImage<unsigned char> &Iresult)
{
    double meanFilter[]{1, 2, 1};
    double derivateFilter[]{-1, 0, 1};
    vpMatrix GradX(3, 3);
    vpMatrix GradY(3, 3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            GradX[i][j] = meanFilter[i] * derivateFilter[j]; // calcul GradX
            GradY[i][j] = meanFilter[j] * derivateFilter[i]; // Calcul gradY
        }
    }

    filtrage2D(I0, IGradX, GradX); // Aplication du filtre à l'image d'origine
    filtrage2D(I0, IGradY, GradY); // Aplication du filtre à l'image d'origine

    for (int i = 0; i < I0.getHeight(); i++)
        for (int j = 0; j < I0.getWidth(); j++)
        {
            Iresult[i][j] = (int)sqrt(pow(IGradX[i][j], 2) + pow(IGradY[i][j], 2)) / 2; // Conversion de l'image obtenue
        }
}

// Rotation du gabarit de k*PI/4
vpImage<unsigned char> rotation_masque(const vpImage<unsigned char> &gabarit, int k)
{

    vpImage<unsigned char> rotated(gabarit.getHeight(), gabarit.getWidth(), 0);

    int half_h = gabarit.getHeight() / 2;
    int half_w = gabarit.getWidth() / 2;
    int offseti, offsetj;

    for (int i = -half_h; i <= half_h; i++)
    {
        for (int j = -half_h; j <= half_h; j++)
        {
            offseti = round(i * cos(k * M_PI / 4) - j * sin(k * M_PI / 4));
            offsetj = round(i * sin(k * M_PI / 4) + j * cos(k * M_PI / 4));
            rotated[offseti + half_h][offsetj + half_h] = gabarit[i + half_h][j + half_h];
        }
    }
    return rotated;
}

void GradientEtMaxima()
{
    // Image originale
    vpImage<unsigned char> I0;
    vpImageIo::read(I0, "../images/contours/elephants.pgm");
    afficheImage(I0, 100, 100, "Image originale");
    const int h = I0.getHeight();
    const int w = I0.getWidth();

    // Tableau des 8 directions
    double angles[8];
    for (int i = 0; i < 8; i++)
    {
        angles[i] = i * M_PI / 4;
    }
    // element structurant non nul à l'horizontale: on peut lui appliquer une roation avec rotation_masque
    vpImage<unsigned char> L(3, 3, 1);
    L[0][0] = 0;
    L[0][1] = 0;
    L[0][2] = 0;
    L[1][0] = 1;
    L[1][1] = 0;
    L[1][2] = 1;
    L[2][0] = 0;
    L[2][1] = 0;
    L[2][2] = 0;

    // Stockage du gradient
    vpImage<double> IgradX(h, w, 0);
    vpImage<double> IgradY(h, w, 0);
    vpImage<unsigned char> Iresult(h, w, 0);
    gradient(I0, IgradX, IgradY, Iresult);
    afficheImage(Iresult, 100, 100, "Gradient");

    for (int i = 0; i < I0.getHeight(); i++)
    {
        for (int j = 0; j < I0.getWidth(); j++)
        {
            double imax = 0;
            double resultMax = 0;
            for (int k = 0; k < 8; k++)
            {
                double angleTempo = k * M_PI / 4;
                double resultTemp = IgradX[i][j] * std::cos(angleTempo) + IgradY[i][j] * std::sin(angleTempo);
                if (resultTemp < 0)
                {
                    resultTemp += 2 * M_PI;
                }
                if (resultTemp > resultMax)
                {
                    resultMax = resultTemp;
                    imax = k;
                }
            }
            vpImage<unsigned char> filtreTempo = rotation_masque(L, imax);
            for (int Ei = 0; Ei < filtreTempo.getHeight(); Ei++)
            {
                for (int Ej = 0; Ej < filtreTempo.getWidth(); Ej++)
                {
                    int itempo = i - filtreTempo.getHeight() / 2 + Ei;
                    int jtempo = j - filtreTempo.getWidth() / 2 + Ej;
                    if (itempo >= 0 && itempo < I0.getHeight() && jtempo >= 0 && jtempo <= I0.getWidth())
                    {
                        if (filtreTempo[Ei][Ej] == 1)
                        {
                            if (Iresult[i][j] < Iresult[itempo][jtempo])
                            {
                                Iresult[i][j] = 0;
                                std::cout << "prout" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    afficheImage(Iresult, 100, 100, "algo");

}

void LaplacienPassageZero(int seuil)
{
    vpImage<unsigned char> L(3, 3, 1);
    L[0][0] = 0;
    L[0][1] = 0;
    L[0][2] = 0;
    L[1][0] = 1;
    L[1][1] = 0;
    L[1][2] = 1;
    L[2][0] = 0;
    L[2][1] = 0;
    L[2][2] = 0;

    // Laplacien
    vpMatrix FIltreLaplacien(3, 3, -0.125); // Laplacien 8-voisins
    FIltreLaplacien[1][1] = 1;

    // Image originale
    vpImage<unsigned char> I0;
    vpImageIo::read(I0, "../images/contours/elephants.pgm");
    afficheImage(I0, 100, 100, "Image originale");
    const int h = I0.getHeight();
    const int w = I0.getWidth();

    // Calcul du laplacien sauvé dans Iresult
    vpImage<double> Itemp(I0.getHeight(), I0.getWidth());
    vpImage<unsigned char> Iresultlaplacien(I0.getHeight(), I0.getWidth());
    filtrage2D(I0, Itemp, FIltreLaplacien);
    // Pour voir le laplacien
    for (int i = 0; i < I0.getHeight(); i++)
    {
        for (int j = 0; j < I0.getWidth(); j++)
        {
            Iresultlaplacien[i][j] = (unsigned char)(Itemp[i][j] + 128);
        }
    }
    afficheImage(Iresultlaplacien, 100, 100, "Image lapalcien");

    vpImage<unsigned char> IresultAlgo(I0.getHeight(), I0.getWidth());
    for (int i = 0; i < I0.getHeight(); i++)
    {
        for (int j = 0; j < I0.getWidth(); j++)
        {
            for (int k = 0; k < 4; k++)
            {
                vpImage<unsigned char> filtreTempo = rotation_masque(L, k);
                double mini = 0;
                double maxi = 0;
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        if (i + l - 1 >= 0 && i + l - 1 < I0.getHeight() && j + m - 1 >= 0 && j + m - 1 < I0.getWidth())
                        {
                            if (filtreTempo[l][m] == 1 && Itemp[i + l - 1][j + m - 1] < 0)
                            {
                                mini = Itemp[i + l - 1][j + m - 1];
                            }
                            if (filtreTempo[l][m] == 1 && Itemp[i + l - 1][j + m - 1] > 0)
                            {
                                maxi = Itemp[i + l - 1][j + m - 1];
                            }
                        }
                        if (maxi * mini < 0 && maxi - mini > seuil)
                        {
                            IresultAlgo[i][j] = 255;
                        }
                    }
                }
            }
        }
    }

    afficheImage(IresultAlgo, 100, 100, "Image algo");

}

void DOG()
{
    // Image originale
    vpImage<unsigned char> I0;
    vpImageIo::read(I0, "../images/contours/elephants.pgm");
    afficheImage(I0, 100, 100, "Image originale");
    const int h = I0.getHeight();
    const int w = I0.getWidth();
    // def des sigma
    double sigma1 = 1.6;
    double sigma2 = 1.0;
    int n = 7;
    int moitie = n / 2;

    // créer les gaussiens
    vpMatrix Gaussien1(n, n);
    vpMatrix Gaussien2(n, n);
    vpMatrix GaussienSous(n, n);
    for (int i = -moitie; i < moitie+1; i++)
    {
        for (int j = -moitie; j < moitie+1; j++)
        {
            Gaussien1[i + moitie][j + moitie] = (1 / (2 * M_PI * std::pow(sigma1, 2)) )* std::exp(-(i * i + j * j) / (2 * sigma1 * sigma1));
            Gaussien2[i + moitie][j + moitie] = 1 / (2 * M_PI * std::pow(sigma2, 2)) * std::exp(-(i * i + j * j) / (2 * sigma2 * sigma2));
            GaussienSous[i + moitie][j + moitie] = Gaussien2[i + moitie][j + moitie] -Gaussien1[i + moitie][j + moitie] ;
        }
    }
    std::cout<<Gaussien1<<std::endl;
    
    vpImage<double> Itemp(h, w);
    vpImage<unsigned char> Iresult(h, w);
    filtrage2D(I0, Itemp, GaussienSous);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            Iresult[i][j]=(unsigned char)Itemp[i][j];
        }
    }
    
    afficheImage(Iresult, 100, 100, "Image algo");
            
}

//
void evaluation(){
      // Image filtré
    vpImage<unsigned char> I0;
    vpImageIo::read(I0, "../images/contours/elephantsGrad.pgm");
    afficheImage(I0, 100, 100, "Image originale");
    const int h = I0.getHeight();
    const int w = I0.getWidth();

    // Image reference
    vpImage<unsigned char> Iref;
    vpImageIo::read(Iref, "../images/contours/test.pgm");
    afficheImage(Iref, 100, 100, "Image originale");

    //Valeurs utiles
    double contours_detected=0;
    double contours_reference=0;
    double contours_corrects=0;
    double faux_positifs=0;
    double faux_negatifs=0;

    for(int i=0;i<h;i++){
        for(int j=0;j<w;j++){
            if(I0[i][j]!=0){
                contours_detected+=1;
            }
            if(Iref[i][j]==0){
                contours_reference+=1;
            }
            if(Iref[i][j]==0 && I0[i][j]!=0){
                contours_corrects+=1;
            }
            if(I0[i][j]!=0 && Iref[i][j]==255 ){
                faux_positifs+=1;
            }
            if(I0[i][j]!=0 && Iref[i][j]==0 ){
                faux_negatifs+=1;
            }
        }
    }
    double performance=contours_corrects/(contours_corrects+faux_positifs+faux_negatifs);
    double taux_de_positif=faux_positifs/(contours_corrects+faux_positifs+faux_negatifs);
    double taux_faux_negatifs=faux_negatifs/(contours_corrects+faux_positifs+faux_negatifs);
    std::cout<<"performance: "<<performance<<std::endl;
    std::cout<<"taux de positif: "<<taux_de_positif<<std::endl;
    std::cout<<"taux faux negatifs: "<<taux_faux_negatifs<<std::endl;

}

int main(int argc, char **argv)
{

    cout << "TIA TP : DETECTION DE CONTOURS " << endl;
    cout << "--" << endl;

    // GradientEtMaxima();
   // LaplacienPassageZero(25);
    //DOG();
    evaluation();

    cout << "Fin du programme " << endl;
    return (0);
}
