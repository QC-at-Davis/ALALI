from distutils.core import setup
setup(
  name = 'ALALI',
  packages = ['ALALI'],
  version = '1.0.0',
  license='gpl-3.0',
  description = 'Automated Loading assuming Atomic Level Interactions (ALALI): Mapping molecules to quantum computers',
  author = 'Jack Goon, John Long, Samarth Sandeep',
  author_email = 'jackbgoon@gmail.com',
  url = 'https://github.com/QC-at-Davis/ALALI',
  download_url = 'https://github.com/QC-at-Davis/ALALI/archive/1.0.0.tar.gz',
  keywords = ['CHEMISTRY', 'QUANTUM', 'QUANTUM CHEMISTRY', 'QUANTUM COMPUTING', 'ATOM', 'QUBIT', 'MAP', 'GRAPH'],   # Keywords that define your package best
  install_requires=[            # pip dependencies, not sure what to do w RDKit
          'networkx',
          'pytket',
          'pytket-qiskit' #rdkit required, must be installed as a prerequisite
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',   
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
