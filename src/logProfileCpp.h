Eigen::MatrixXd logProfileCpp(const Eigen::VectorXd theta0,
                              const int CovStructure,
                              const Eigen::MatrixXd loctim,
                              const Eigen::MatrixXd DS,
                              const Eigen::MatrixXd DT,
                              const Eigen::MatrixXd Xprime,
                              const Eigen::VectorXd Zprime,
                              const Eigen::VectorXd truncate,
                              const double lambda,
                              const int pterm,
                              const Eigen::VectorXd DDnew
                              );
