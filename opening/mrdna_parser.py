import argparse

def make_parser():
    import argparse
    #specific options for staplehedra
    parser = argparse.ArgumentParser(description = "Make poly")


    parser.add_argument('-o','--output-prefix', type=str, default=None,
                        help="Name for your job's output")
    parser.add_argument('-d','--directory', type=str,  default=None,
                        help='Directory for simulation; does not need to exist yet')
    parser.add_argument('-g','--gpu', type=int, default=0,
                        help='GPU for simulation; check nvidia-smi for availability')
    parser.add_argument('--debye-length', type=float, default=None,
                        help='Adjust the electrostatic repulsion matching osmotic pressure data from Rau and Parsegian under 25 mM MgCl2 condition with a Debye-Hueckel correction from the default 11.1 Angstroms.')
    parser.add_argument('--temperature', type=float, default=295,
                        help='Temperature in Kelvin.')
    parser.add_argument('--sequence-file', type=str, default=None,
                        help='Sequence of longest strand.')
    parser.add_argument('--output-period', type=float, default=1e4,
                        help='Simulation steps between DCD frames')
    # parser.add_argument('--minimization-steps', type=float, default=0,
    #                     help='Simulation steps between DCD frames')
    parser.add_argument('--coarse-local-twist', action='store_true',
                        help='Use local twist for coarsest representation?')
    parser.add_argument('--fix-linking-number', action='store_true',
                        help='Fix the linking number so that it cannot change once local twist is introduced?')
    parser.add_argument('--coarse-steps', type=float, default=1e7,
                        help='Simulation steps for coarse model (200 fs/step)')
    parser.add_argument('--fine-steps', type=float, default=1e7,
                        help='Simulation steps for fine model (50 fs/step)')

    parser.add_argument('--oxdna-steps', type=float, default=None,
                        help='Perform an oxDNA simulation instead of creating atomic model')
    parser.add_argument('--oxdna-output-period', type=float, default=1e4,
                        help='Simulation steps between oxDNA configuration and energy output')

    parser.add_argument('--backbone-scale', type=float, default=1.0,
                        help='Factor to scale DNA backbone in atomic model; try 0.25 to avoid clashes for atomistic simulations')
    parser.add_argument('--debug', action='store_true',
                        help='Run through the python debugger?')

    parser.add_argument('--draw-cylinders', action='store_true',
                        help='Whether or not to draw the cylinders')
    parser.add_argument('--draw-tubes', action='store_true',
                        help='Whether or not to draw the tubes')

    return parser