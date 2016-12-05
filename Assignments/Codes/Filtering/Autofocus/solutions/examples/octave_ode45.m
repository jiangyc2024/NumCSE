%# Copyright (C) 2006-2012, Thomas Treichl <treichl@users.sourceforge.net>
%# Copyright (C) 2013, Roberto Porcu' <roberto.porcu@polimi.it>
%# OdePkg - A package for solving ordinary differential equations and more
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%#
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with this program; If not, see <http://www.gnu.org/licenses/>.

%# -*- texinfo -*-
%# @deftypefn {Command} {[@var{sol}] =} ode45 (@var{@@fun}, @var{slot}, @var{init}, [@var{opt}], [@var{par1}, @var{par2}, @dots{}])
%# @deftypefnx {Command} {[@var{t}, @var{y}, [@var{xe}, @var{ye}, @var{ie}]] =} ode45 (@var{@@fun}, @var{slot}, @var{init}, [@var{opt}], [@var{par1}, @var{par2}, @dots{}])
%# 
%# This function file can be used to solve a set of non--stiff ordinary differential equations (non--stiff ODEs) with the well known explicit Dormand-Prince method of order 4.
%# 
%# This function can be called with two output arguments: @var{t} and @var{y}. Variable @var{t} is a column vector and contains the time stamps, instead @var{y} is a matrix in which each column refers to a different unknown of the problem and the rows number is the same of @var{t} rows number so that each row of @var{y} contains the values of all unknowns at the time value contained in the corresponding row in @var{t}.
%# 
%# The first input argument must be a function_handle or an inline function that defines the set of ODE: @code{y' = f(t,y)}. As described above, this function must take two input arguments, where the first is the time and the second the unknowns, and must have just one output argument.
%# 
%# The second input argument must contain time informations. Usually it should be a vector with at least two elements which define the initial and the final time instants; if the elements are more than two, then the solution will be evaluated also at these intermediate time instants unless the integrate function called is the @command{integrate_n_steps}. If there is only one time value, then it will give an error unless the options structure has no empty fields named @var{"TimeStepNumber"} and @var{"TimeStepSize"}. If the option @var{"TimeStepSize"} is not empty, then the stepper called will be @command{integrate_const}, if also @var{"TimeStepNumber"} is not empty it will be called the integrate function @command{integrate_n_steps}, otherwise it will be called @command{integrate_adaptive}. For this last possibility the user can set the tolerance for the timestep computation by setting a value to the option @var{"Tau"}, that as default value has @math{1.e-6}.
%# 
%# The third input argument must contain the initial value for the unknown. If this is a vector then the solution @var{y} will be a matrix in which each column is the solution for the corresponding initial value in @var{init}.
%# 
%# The fourth input argument is not mandatory and it should contain a structure with valid ODE fields.
%# 
%# For example, solve an anonymous implementation of the Van der Pol equation
%# @example
%# fvdp = @@(t,y) [y(2); (1 - y(1)^2) * y(2) - y(1)];
%# [T,Y] = ode45 (fvdp, [0 20], [2 0]);
%# @end example
%# @end deftypefn
%#
%# @seealso{odepkg}

%# ChangeLog:
%#   20010703 the function file "ode45.m" was written by Marc Compere
%#     under the GPL for the use with this software. This function has been
%#     taken as a base for the following implementation.
%#   20060810, Thomas Treichl
%#     This function was adapted to the new syntax that is used by the
%#     new OdePkg for Octave and is compatible to Matlab's ode45.

function [varargout] = ode45 (vfun, vslot, vinit, vodeoptions, varargin)

  vorder = 4;
  vsolver = 'ode45';

  if (nargin == 0) %# Check number and types of all input arguments
    help (vsolver);
    error ('OdePkg:InvalidArgument', ...
      'number of input arguments must be greater than zero');
  end

  if (nargin < 3)
    print_usage;
  elseif (nargin < 4)
    vodeoptions = odeset;
  elseif (nargin > 4)
    warning('OdePkg:InvalidArgument', ...
      'all inputs after fourth input argument will be ignored');
  end

  if (~isvector (vslot) || ~isnumeric(vslot))
    error ('OdePkg:InvalidArgument', ...
      'second input argument must be a valid vector');
  end

  if(length(vslot)<2 && (isempty(vodeoptions.TimeStepSize) || isempty(vodeoptions.TimeStepNumber) ))
    error ('OdePkg:InvalidArgument', ...
      'second input argument must be a valid vector');
  end

  if (~isvector (vinit) || ~isnumeric (vinit))
    error ('OdePkg:InvalidArgument', ...
      'third input argument must be a valid numerical value');
  end
  vinit = vinit(:);

  if ~(isa (vfun, 'function_handle') || isa (vfun, 'inline'))
    error ('OdePkg:InvalidArgument', ...
      'first input argument must be a valid function handle');
  end

  if(~isempty(vodeoptions))
    if(isstruct(vodeoptions))
      warning('OdePkg:InvalidArgument', ...
        'fourth input argument should be a valid ODE structure. It will be ignored');
      vodeoptions = odeset;
    else
      try
        odepkg_structure_check(vodeoptions,vsolver);
      catch
        vodeoptions = odeset();
        warning('OdePkg:InvalidArgument', ...
          'fourth input argument should be a valid ODE structure. It will be ignored');
      end
    end
  end

  %# Start preprocessing, have a look which options are set in
  %# vodeoptions, check if an invalid or unused option is set
  if(isempty(vodeoptions.TimeStepNumber) && isempty(vodeoptions.TimeStepSize))
    integrate_func = 'adaptive';
  elseif(~isempty(vodeoptions.TimeStepNumber) && ~isempty(vodeoptions.TimeStepSize))
    integrate_func = 'n_steps';
  elseif(isempty(vodeoptions.TimeStepNumber) && ~isempty(vodeoptions.TimeStepSize))
    integrate_func = 'const';
  else
    warning ('OdePkg:InvalidArgument', ...
      'assuming an adaptive integrate function');
    integrate_func = 'adaptive';
  end  

  %# Get the default options that can be set with 'odeset' temporarily
  vodetemp = odeset;

  %# Implementation of the option RelTol has been finished. This option
  %# can be set by the user to another value than default value.
  if (~isequal (vodeoptions.RelTol, vodetemp.RelTol))
    warning ('OdePkg:InvalidArgument', ...
      'option "RelTol" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option AbsTol has been finished. This option
  %# can be set by the user to another value than default value.
  if (~isequal (vodeoptions.AbsTol, vodetemp.AbsTol))
    warning ('OdePkg:InvalidArgument', ...
      'option "AbsTol" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option NormControl has been finished. This
  %# option can be set by the user to another value than default value.
  if (~isequal (vodeoptions.NormControl, vodetemp.NormControl))
    warning ('OdePkg:InvalidArgument', ...
      'option "NormControl" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option NonNegative has been finished. This
  %# option can be set by the user to another value than default value.
  if (~isequal (vodeoptions.NonNegative, vodetemp.NonNegative))
    warning ('OdePkg:InvalidArgument', ...
      'option "NonNegative" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option OutputFcn has been finished. This
  %# option can be set by the user to another value than default value.
  if (~isequal (vodeoptions.OutputFcn, vodetemp.OutputFcn))
    warning ('OdePkg:InvalidArgument', ...
      'option "OutputFcn" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option OutputSel has been finished. This
  %# option can be set by the user to another value than default value.
  if (~isequal (vodeoptions.OutputSel, vodetemp.OutputSel))
    warning ('OdePkg:InvalidArgument', ...
      'option "OutputSel" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option OutputSave has been finished. This
  %# option can be set by the user to another value than default value.
  if (~isequal (vodeoptions.OutputSave, vodetemp.OutputSave))
    warning ('OdePkg:InvalidArgument', ...
      'option "OutputSave" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option Refine has been finished. This option
  %# can be set by the user to another value than default value.
  if (~isequal (vodeoptions.Refine, vodetemp.Refine))
    warning ('OdePkg:InvalidArgument', ...
      'option "Refine" will be ignored by this solver');
  end

  %# Implementation of the option Stats has been finished. This option
  %# can be set by the user to another value than default value.

  %# Implementation of the option InitialStep has been finished. This
  %# option can be set by the user to another value than default value.
  if (isempty (vodeoptions.InitialStep) && strcmp(integrate_func,'adaptive'))
    vodeoptions.InitialStep = starting_stepsize(vorder,vfun,vslot(1),vinit);
    warning ('OdePkg:InvalidArgument', ...
      'option "InitialStep" not set, estimated value %f is used', vodeoptions.InitialStep);
  elseif(isempty (vodeoptions.InitialStep))
    vodeoptions.InitialStep = odeget(vodeoptions,'TimeStepSize');
  end

  %# Implementation of the option MaxStep has been finished. This option
  %# can be set by the user to another value than default value.
  if (~isequal (vodeoptions.MaxStep, vodetemp.MaxStep))
    warning ('OdePkg:InvalidArgument', ...
      'option "MaxStep" will be ignored, not yet implemented functionality');
  end

  %# Implementation of the option Events has been finished. This option
  %# can be set by the user to another value than default value.
  if (~isequal (vodeoptions.Events, vodetemp.Events))
    warning ('OdePkg:InvalidArgument', ...
      'option "Events" will be ignored by this solver');
  end

  %# The options 'Jacobian', 'JPattern' and 'Vectorized' will be ignored
  %# by this solver because this solver uses an explicit Runge-Kutta
  %# method and therefore no Jacobian calculation is necessary
  if (~isequal (vodeoptions.Jacobian, vodetemp.Jacobian))
    warning ('OdePkg:InvalidArgument', ...
      'option "Jacobian" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.JPattern, vodetemp.JPattern))
    warning ('OdePkg:InvalidArgument', ...
      'option "JPattern" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.Vectorized, vodetemp.Vectorized))
    warning ('OdePkg:InvalidArgument', ...
      'option "Vectorized" will be ignored, not yet implemented functionality');
  end

  if (~isequal (vodeoptions.NewtonTol, vodetemp.NewtonTol))
    warning ('OdePkg:InvalidArgument', ...
      'option "NewtonTol" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.MaxNewtonIterations,vodetemp.MaxNewtonIterations))
    warning ('OdePkg:InvalidArgument', ...
      'option "MaxNewtonIterations" will be ignored by this solver');
  end

  %# Implementation of the option Mass has been finished. This option
  %# can be set by the user to another value than default value.
  if (~isequal (vodeoptions.Mass, vodetemp.Mass))
    warning ('OdePkg:InvalidArgument', ...
      'option "Mass" will be ignored by this solver');
  end

  %# Implementation of the option MStateDependence has been finished.
  %# This option can be set by the user to another value than default
  %# value. 
  if (~isequal (vodeoptions.MStateDependence, vodetemp.MStateDependence))
    warning ('OdePkg:InvalidArgument', ...
      'option "MStateDependence" will be ignored by this solver');
  end

  %# Other options that are not used by this solver. Print a warning
  %# message to tell the user that the option(s) is/are ignored.
  if (~isequal (vodeoptions.MvPattern, vodetemp.MvPattern))
    warning ('OdePkg:InvalidArgument', ...
      'option "MvPattern" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.MassSingular, vodetemp.MassSingular))
    warning ('OdePkg:InvalidArgument', ...
      'option "MassSingular" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.InitialSlope, vodetemp.InitialSlope))
    warning ('OdePkg:InvalidArgument', ...
      'option "InitialSlope" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.MaxOrder, vodetemp.MaxOrder))
    warning ('OdePkg:InvalidArgument', ...
      'option "MaxOrder" will be ignored by this solver');
  end

  if (~isequal (vodeoptions.BDF, vodetemp.BDF))
    warning ('OdePkg:InvalidArgument', ...
      'option "BDF" will be ignored by this solver');
  end

  %# Starting the initialisation of the core solver ode45 

  SubOpts = vodeoptions;

  SobOpts.Tau = odeget(vodeoptions,'Tau',1.e-6,'fast_not_empty');

  switch integrate_func
    case 'adaptive'
      SobOpts.Tau = odeget(vodeoptions,'Tau',1.e-6,'fast_not_empty');
      [t,x] = integrate_adaptive(@runge_kutta_45,vorder,vfun,vslot,vinit,SubOpts);
    case 'n_steps'
      [t,x] = integrate_n_steps(@runge_kutta_45,vfun,vslot(1),vinit,vodeoptions.TimeStepSize,vodeoptions.TimeStepNumber,SubOpts);
    case 'const'
      [t,x] = integrate_const(@runge_kutta_45,vfun,vslot,vinit,vodeoptions.TimeStepSize,SubOpts);
  end

  if(nargout==1)
    varargout{1}.x = t;
    varargout{1}.y = x;
    varargout{1}.solver = solver;
  elseif(nargout==2)
    varargout{1} = t;
    varargout{2} = x;
  elseif(nargout==5)
    varargout{1} = t;
    varargout{2} = x;
    varargout{3} = [];
    varargout{4} = [];
    varargout{5} = [];
  end

end

%! # We are using the "Van der Pol" implementation for all tests that
%! # are done for this function.
%! # For further tests we also define a reference solution (computed at high accuracy)
%!function [ydot] = fpol (vt, vy) %# The Van der Pol
%!  ydot = [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!endfunction
%!function [vref] = fref ()         %# The computed reference sol
%!  vref = [0.32331666704577, -1.83297456798624];
%!endfunction
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!error %# ouput argument
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  B = ode45 (1, [0 25], [3 15 1]);
%!error %# input argument number one
%!  [vt, vy] = ode45 (1, [0 25], [3 15 1]);
%!error %# input argument number two
%!  [vt, vy] = ode45 (@fpol, 1, [3 15 1]);
%!test %# two output arguments
%!  [vt, vy] = ode45 (@fpol, [0 2], [2 0]);
%!  assert ([vt(end), vy(end,:)], [2, fref], 1e-2);
%!test %# anonymous function instead of real function
%!  fvdb = @(vt,vy) [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!  [vt, vy] = ode45 (fvdb, [0 2], [2 0]);
%!  assert ([vt(end), vy(end,:)], [2, fref], 1e-2);
%!test %# extra input arguments passed through
%!  [vt, vy] = ode45 (@fpol, [0 2], [2 0], 12, 13, 'KL');
%!  assert ([vt(end), vy(end,:)], [2, fref], 1e-2);
%!test %# empty OdePkg structure *but* extra input arguments
%!  vopt = odeset;
%!  [vt, vy] = ode45 (@fpol, [0 2], [2 0], vopt, 12, 13, 'KL');
%!  assert ([vt(end), vy(end,:)], [2, fref], 1e-2);
%!test %# strange OdePkg structure
%!  vopt = struct ('foo', 1);
%!  [vt, vy] = ode45 (@fpol, [0 2], [2 0], vopt);
%!test %# Solve vdp in fixed step sizes
%!  vopt = odeset('TimeStepSize',0.1);
%!  [vt, vy] = ode45 (@fpol, [0,2], [2 0], vopt);
%!  assert (vt(:), [0:0.1:2]', 1e-2);
%!test %# Solve another anonymous function below zero
%!  vref = [0, 14.77810590694212];
%!  [vt, vy] = ode45 (@(t,y) y, [-2 0], 2);
%!  assert ([vt(end), vy(end,:)], vref, 1e-1);
%!test %# InitialStep option
%!  vopt = odeset ('InitialStep', 1e-8);
%!  [vt, vy] = ode45 (@fpol, [0 0.2], [2 0], vopt);
%!  assert ([vt(2)-vt(1)], [1e-8], 1e-9);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
