# Modified the controller and launcher in order to have r command be executed from the singularity container

async_launcher_class <- R6::R6Class(
  classname = "custom_launcher_class",
  inherit = crew::crew_class_launcher,
  public = list(
    launch_worker = function(call, name, launcher, worker, instance) {
      self$async$eval(
        command = list(pid = process$new(bin, args = c("-e", call))$get_pid()),
        data = list(bin = "singularity exec /mnt/beegfs/idevillasante/apps/rocker/images/meth-dev.sif r -e ", call = call),
        packages = "processx"
      )
    },
    terminate_worker = function(handle) {
      self$async$eval(
        command = ps::ps_kill(p = ps::ps_handle(pid = pid)),
        data = list(pid = handle$data$pid)
      )
    }
  )
)

crew_controller_async <- function(
    name = "async controller name",
    workers = 1L,
    host = "127.0.0.1",
    port = NULL,
    tls = crew::crew_tls(mode = "none"),
    seconds_interval = 0.5,
    seconds_timeout = 5,
    seconds_launch = 30,
    seconds_idle = Inf,
    seconds_wall = Inf,
    tasks_max = Inf,
    tasks_timers = 0L,
    reset_globals = TRUE,
    reset_packages = FALSE,
    reset_options = FALSE,
    garbage_collection = FALSE,
    launch_max = 5L,
    processes = NULL # Number of local async daemons for worker launches etc.
) {
  client <- crew::crew_client(
    name = name,
    workers = workers,
    host = host,
    port = port,
    tls = tls,
    seconds_interval = seconds_interval,
    seconds_timeout = seconds_timeout
  )
  launcher <- async_launcher_class$new(
    name = name,
    seconds_interval = seconds_interval,
    seconds_launch = seconds_launch,
    seconds_idle = seconds_idle,
    seconds_wall = seconds_wall,
    tasks_max = tasks_max,
    tasks_timers = tasks_timers,
    reset_globals = reset_globals,
    reset_packages = reset_packages,
    reset_options = reset_options,
    garbage_collection = garbage_collection,
    launch_max = launch_max,
    tls = tls,
    processes = processes
  )
  controller <- crew::crew_controller(
    client = client,
    launcher = launcher
  )
  controller$validate()
  controller
}

