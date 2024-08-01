float PiRange(float deltaPhi);
bool IsClusterClose(float projPhi, float projEta, float clusterPhi, float clusterEta, float minR);
bool HasOverlap(float projPhi, float projEta, std::vector<float> towerPhi, std::vector<float> towerEta, float minDelta);
void DrawEllipseWithTGraph(double x_center, double y_center, double radius, int color);
void DrawVLineWithTGraph(double x_center, double y_min, double y_max, int color);
void DrawHLineWithTGraph(double y_center, double x_min, double x_max, int color);

int _EMCAL_NETA = 96;
int _EMCAL_NPHI = 256;
int _HCAL_NETA = 24;
int _HCAL_NPHI = 64;

float min_EMCal_E = 1.5;
float min_HCal_E = 0.2;

float emcal_radius = 100.70;//(1-(-0.077))*93.5
float hcal_radius = 177.423;

float PiRange(float deltaPhi)
{
  if(deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  if(deltaPhi < -M_PI) deltaPhi += 2*M_PI;

  return deltaPhi;
}
bool IsClusterClose(float projPhi, float projEta, float clusterPhi, float clusterEta, float minR)
{
  if(isnan(projPhi) || isnan(projEta))
  {
    return false;
  }
  float deltaEta = projEta - clusterEta;
  float deltaPhi = projPhi - clusterPhi;
  deltaPhi = PiRange(deltaPhi);
  float deltaR = sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

  if(deltaR < minR)
  {
    //std::cout << "Cluster deltaPhi: " << deltaPhi << " deltaEta: " << deltaEta << " deltaR: " << deltaR << std::endl;
    return true;
  }

  return false;
}
bool HasOverlap(float projPhi, float projEta, std::vector<float> towerPhi, std::vector<float> towerEta, float minDelta)
{
  for(unsigned int i = 0; i < towerPhi.size(); i++)
  {
    float deltaEta = std::fabs(projEta - towerEta.at(i));
    float deltaPhi = projPhi - towerPhi.at(i);
    deltaPhi = std::fabs(PiRange(deltaPhi));
    if((deltaPhi < minDelta) && (deltaEta < minDelta))
    {
      return true;
    }
  }

  return false;
}

void DrawEllipseWithTGraph(double x_center, double y_center, double radius, int color) {
    const int n_points = 1000; // Number of points to approximate the ellipse
    TGraph *gr_ellipse = new TGraph(n_points);
    for (int i = 0; i < n_points; ++i) {
        double theta = 2 * TMath::Pi() * i / n_points;
        double x = x_center + radius * TMath::Cos(theta);
        double y = y_center + radius * TMath::Sin(theta);
        gr_ellipse->SetPoint(i, x, y);
    }
    gr_ellipse->SetLineColor(color);
    gr_ellipse->SetLineWidth(2);
    gr_ellipse->Draw("L same");
}

void DrawVLineWithTGraph(double x_center, double y_min, double y_max, int color) {
    const int n_points = 1000; // Number of points to approximate the line
    TGraph *gr_line = new TGraph(n_points);
    for (int i = 0; i < n_points; ++i) {
        double x = x_center;
        double y = y_min + i * (y_max - y_min) / (double) n_points;
        gr_line->SetPoint(i, x, y);
    }
    gr_line->SetLineColor(color);
    gr_line->SetLineWidth(2);
    gr_line->Draw("L same");
}
void DrawHLineWithTGraph(double y_center, double x_min, double x_max, int color) {
    const int n_points = 1000; // Number of points to approximate the line
    TGraph *gr_line = new TGraph(n_points);
    for (int i = 0; i < n_points; ++i) {
        double x = x_min + i * (x_max - x_min) / (double) n_points;
        double y = y_center;
        gr_line->SetPoint(i, x, y);
    }
    gr_line->SetLineColor(color);
    gr_line->SetLineWidth(2);
    gr_line->Draw("L same");
}
