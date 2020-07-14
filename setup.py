from distutils.core import setup
setup(
  name = 'ALALI',
  packages = ['ALALI'],
  version = '0.0.1',
  license='gpl-3.0',
  description = 'Automated Loading assuming Atomic Level Interactions (ALALI): Mapping molecules to quantum computers',
  author = 'Jack Goon, John Long, Samarth Sandeep',
  author_email = 'jackbgoon@gmail.com',
  url = 'https://github.com/QC-at-Davis/ALALI',
  download_url = 'https://github.com/QC-at-Davis/ALALI/archive/0.0.1.tar.gz', #Need to upload to github before getting download link
  keywords = ['CHEMISTRY', 'QUANTUM', 'QUANTUM CHEMISTRY', 'QUANTUM COMPUTING', 'ATOM', 'QUBIT', 'MAP', 'GRAPH'],   # Keywords that define your package best
  install_requires=[
          'networkx',
          'pytket',
          'pytket-qiskit',
          'matplotlib',
          'rdkit @ https://github.com/rdkit/rdkit/archive/master.zip#egg=rdkit' #not sure on this
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU General Public License (GPL)',   
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
  ],
)
