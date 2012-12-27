// Dear emacs, this is -*-c++-*-
// 
// Andrei Gaponenko, 2009

#ifndef GETHISTO_H
#define GETHISTO_H

#include <string>

class TH1;

TH1* getHisto(const std::string& filename, const std::string& histoname);

#endif/*GETHISTO_H*/

