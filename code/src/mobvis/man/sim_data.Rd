% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_functions.R
\name{sim_get_region}
\alias{sim_get_region}
\alias{sim_get_raster}
\alias{sim_get_cellplan}
\alias{sim_get_signal_strength}
\alias{sim_get_trajectory_data}
\alias{sim_get_trajectory_data_TA}
\alias{sim_get_trajectory_routes}
\alias{sim_devices_at_t}
\alias{sim_get_prob}
\title{Get simulator data}
\usage{
sim_get_region(sim)

sim_get_raster(sim)

sim_get_cellplan(sim)

sim_get_signal_strength(sim, rst, cp)

sim_get_trajectory_data(sim, device = NULL)

sim_get_trajectory_data_TA(sim, device = NULL)

sim_get_trajectory_routes(sim, device = NULL)

sim_devices_at_t(sim, t)

sim_get_prob(sim, device, rootfilename, th = 1e-06)
}
\arguments{
\item{sim}{list containing the following elements: \code{simConfig_dir},
\code{resources_dir}, \code{networkEvents_dir}, \code{postLocProb_dir},
\code{mno}, \code{crs}.}

\item{rst}{rst}

\item{cp}{cp}

\item{device}{device id}

\item{t}{time id}

\item{th}{probability threshold. Only probabilities of at least `th` are returned}
}
\description{
Get simulator data
}
