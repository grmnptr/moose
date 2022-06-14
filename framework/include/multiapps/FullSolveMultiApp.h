//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MultiApp.h"

class Executioner;

/**
 * This type of MultiApp will do a full solve when it is asked to take a step.
 */
class FullSolveMultiApp : public MultiApp
{
public:
  static InputParameters validParams();

  FullSolveMultiApp(const InputParameters & parameters);

  virtual void initialSetup() override;

  /// Resets the executioners of the owned applications. Used in cases when
  /// we don't want to reset executioners without reinitializing the whole system
  /// @param needs_reset a vector of flags showing which executioners needs to be reset
  virtual void resetExecutioners(const std::vector<bool> & needs_reset = std::vector<bool>());

  /// Reinitialize applications. This means that the old applications are destroyed and
  /// new ones are created to replace the old.
  /// @param needs_reset a vector of flags showing which executioners needs to be reset
  virtual void reinitApps(const std::vector<bool> & needs_reset = std::vector<bool>());

  virtual bool solveStep(Real dt, Real target_time, bool auto_advance = true) override;

  virtual void finalize() override
  {
    // executioner output on final has been called and we do not need to call it again
  }
  virtual void postExecute() override
  {
    // executioner postExecute has been called and we do not need to call it again
  }

  virtual void backup() override;
  virtual void restore(bool force = true) override;

  virtual bool ignoreDiverge() { return _ignore_diverge; }

  std::vector<Executioner *> & setExecutioners() { return _executioners; }
  const std::vector<Executioner *> & getExecutioners() const { return _executioners; }

private:
  /// Switch to tell executioner to keep going despite app solve not converging
  const bool _ignore_diverge;

  std::vector<Executioner *> _executioners;
};
