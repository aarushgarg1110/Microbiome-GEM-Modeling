FROM python:3.10-slim

# Set environment for CPLEX
ENV CPLEX_HOME=/opt/ibm/cplex
ENV PATH="${CPLEX_HOME}/bin/x86-64_linux:$PATH"
ENV LD_LIBRARY_PATH="${CPLEX_HOME}/lib/x86-64_linux:$LD_LIBRARY_PATH"

# Minimal runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libglib2.0-0 libgomp1 libquadmath0 && \
    rm -rf /var/lib/apt/lists/*

# Copy installed CPLEX files (from host)
COPY cplex /opt/ibm/cplex

# Install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install CPLEX Python bindings (must match Python version)
COPY cplex/python/3.10/x86-64_linux /tmp/cplex_py/
RUN cd /tmp/cplex_py && python setup.py install && rm -rf /tmp/cplex_py

# Copy your modeling script
COPY script.py /app/script.py
WORKDIR /app

CMD ["python", "script.py"]

PS C:\Users\agarg> docker run -it --rm -v "${PWD}:/mnt" debian:bullseye bash
root@bd9f6e7f9f8e:/# exit
exit
PS C:\Users\agarg> cd .\research\Microbiome-GEM-Modeling\
PS C:\Users\agarg\research\Microbiome-GEM-Modeling> docker run -it --rm -v "${PWD}:/mnt" debian:bullseye bash
root@79cf31c19188:/# ls
bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
root@79cf31c19188:/# ls opt
root@79cf31c19188:/# ls mnt
 Dockerfile                 Python_Models2   __pycache__   misc               src
 Matlab_Models              README.md        cplex_stuff   requirements.txt   test_data_input
'ModelBuilder copy.ipynb'   Summary.md       mgPipe        script.py          tests.ipynb
root@79cf31c19188:/# apt update && apt install -y openjdk-11-jre-headless
Get:1 http://deb.debian.org/debian bullseye InRelease [116 kB]
Get:2 http://deb.debian.org/debian-security bullseye-security InRelease [27.2 kB]
Get:3 http://deb.debian.org/debian bullseye-updates InRelease [44.0 kB]
Get:4 http://deb.debian.org/debian bullseye/main amd64 Packages [8066 kB]
Get:5 http://deb.debian.org/debian-security bullseye-security/main amd64 Packages [384 kB]
Get:6 http://deb.debian.org/debian bullseye-updates/main amd64 Packages [18.8 kB]
Fetched 8656 kB in 3s (2963 kB/s)
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
All packages are up to date.
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following additional packages will be installed:
  alsa-topology-conf alsa-ucm-conf ca-certificates ca-certificates-java dbus fontconfig-config fonts-dejavu-core
  java-common libapparmor1 libasound2 libasound2-data libavahi-client3 libavahi-common-data libavahi-common3
  libbrotli1 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libfreetype6 libglib2.0-0 libglib2.0-data libgraphite2-3
  libharfbuzz0b libicu67 libjpeg62-turbo liblcms2-2 libnspr4 libnss3 libpcsclite1 libpng16-16 libsqlite3-0 libxml2
  openssl sensible-utils shared-mime-info ucf xdg-user-dirs
Suggested packages:
  default-dbus-session-bus | dbus-session-bus default-jre libasound2-plugins alsa-utils cups-common liblcms2-utils
  pcscd libnss-mdns fonts-dejavu-extra fonts-ipafont-gothic fonts-ipafont-mincho fonts-wqy-microhei | fonts-wqy-zenhei
  fonts-indic
The following NEW packages will be installed:
  alsa-topology-conf alsa-ucm-conf ca-certificates ca-certificates-java dbus fontconfig-config fonts-dejavu-core
  java-common libapparmor1 libasound2 libasound2-data libavahi-client3 libavahi-common-data libavahi-common3
  libbrotli1 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libfreetype6 libglib2.0-0 libglib2.0-data libgraphite2-3
  libharfbuzz0b libicu67 libjpeg62-turbo liblcms2-2 libnspr4 libnss3 libpcsclite1 libpng16-16 libsqlite3-0 libxml2
  openjdk-11-jre-headless openssl sensible-utils shared-mime-info ucf xdg-user-dirs
0 upgraded, 39 newly installed, 0 to remove and 0 not upgraded.
Need to get 60.6 MB of archives.
After this operation, 254 MB of additional disk space will be used.
Get:1 http://deb.debian.org/debian bullseye/main amd64 libapparmor1 amd64 2.13.6-10 [99.3 kB]
Get:2 http://deb.debian.org/debian bullseye/main amd64 libdbus-1-3 amd64 1.12.28-0+deb11u1 [223 kB]
Get:3 http://deb.debian.org/debian-security bullseye-security/main amd64 libexpat1 amd64 2.2.10-2+deb11u7 [99.2 kB]
Get:4 http://deb.debian.org/debian bullseye/main amd64 dbus amd64 1.12.28-0+deb11u1 [244 kB]
Get:5 http://deb.debian.org/debian bullseye/main amd64 sensible-utils all 0.0.14 [14.8 kB]
Get:6 http://deb.debian.org/debian-security bullseye-security/main amd64 openssl amd64 1.1.1w-0+deb11u3 [859 kB]
Get:7 http://deb.debian.org/debian bullseye/main amd64 ca-certificates all 20210119 [158 kB]
Get:8 http://deb.debian.org/debian-security bullseye-security/main amd64 ucf all 3.0043+deb11u2 [74.3 kB]
Get:9 http://deb.debian.org/debian bullseye/main amd64 alsa-topology-conf all 1.2.4-1 [12.8 kB]
Get:10 http://deb.debian.org/debian bullseye/main amd64 libasound2-data all 1.2.4-1.1 [38.2 kB]
Get:11 http://deb.debian.org/debian bullseye/main amd64 libasound2 amd64 1.2.4-1.1 [356 kB]
Get:12 http://deb.debian.org/debian bullseye/main amd64 alsa-ucm-conf all 1.2.4-2 [28.1 kB]
Get:13 http://deb.debian.org/debian bullseye/main amd64 java-common all 0.72 [14.5 kB]
Get:14 http://deb.debian.org/debian-security bullseye-security/main amd64 libavahi-common-data amd64 0.8-5+deb11u3 [124 kB]
Get:15 http://deb.debian.org/debian-security bullseye-security/main amd64 libavahi-common3 amd64 0.8-5+deb11u3 [59.0 kB]
Get:16 http://deb.debian.org/debian-security bullseye-security/main amd64 libavahi-client3 amd64 0.8-5+deb11u3 [62.7 kB]
Get:17 http://deb.debian.org/debian-security bullseye-security/main amd64 libcups2 amd64 2.3.3op2-3+deb11u9 [351 kB]
Get:18 http://deb.debian.org/debian bullseye/main amd64 libbrotli1 amd64 1.0.9-2+b2 [279 kB]
Get:19 http://deb.debian.org/debian bullseye/main amd64 libpng16-16 amd64 1.6.37-3 [294 kB]
Get:20 http://deb.debian.org/debian-security bullseye-security/main amd64 libfreetype6 amd64 2.10.4+dfsg-1+deb11u2 [418 kB]
Get:21 http://deb.debian.org/debian bullseye/main amd64 fonts-dejavu-core all 2.37-2 [1069 kB]
Get:22 http://deb.debian.org/debian bullseye/main amd64 fontconfig-config all 2.13.1-4.2 [281 kB]
Get:23 http://deb.debian.org/debian bullseye/main amd64 libfontconfig1 amd64 2.13.1-4.2 [347 kB]
Get:24 http://deb.debian.org/debian bullseye/main amd64 liblcms2-2 amd64 2.12~rc1-2 [150 kB]
Get:25 http://deb.debian.org/debian bullseye/main amd64 libjpeg62-turbo amd64 1:2.0.6-4 [151 kB]
Get:26 http://deb.debian.org/debian bullseye/main amd64 libnspr4 amd64 2:4.29-1 [112 kB]
Get:27 http://deb.debian.org/debian-security bullseye-security/main amd64 libsqlite3-0 amd64 3.34.1-3+deb11u1 [797 kB]
Get:28 http://deb.debian.org/debian-security bullseye-security/main amd64 libnss3 amd64 2:3.61-1+deb11u4 [1304 kB]
Get:29 http://deb.debian.org/debian-security bullseye-security/main amd64 libglib2.0-0 amd64 2.66.8-1+deb11u6 [1377 kB]
Get:30 http://deb.debian.org/debian bullseye/main amd64 libgraphite2-3 amd64 1.3.14-1 [81.2 kB]
Get:31 http://deb.debian.org/debian bullseye/main amd64 libharfbuzz0b amd64 2.7.4-1 [1471 kB]
Get:32 http://deb.debian.org/debian bullseye/main amd64 libpcsclite1 amd64 1.9.1-1 [60.2 kB]
Get:33 http://deb.debian.org/debian-security bullseye-security/main amd64 openjdk-11-jre-headless amd64 11.0.27+6-1~deb11u1 [38.3 MB]
Get:34 http://deb.debian.org/debian bullseye/main amd64 ca-certificates-java all 20190909+deb11u1 [15.9 kB]
Get:35 http://deb.debian.org/debian-security bullseye-security/main amd64 libglib2.0-data all 2.66.8-1+deb11u6 [1177 kB]
Get:36 http://deb.debian.org/debian-security bullseye-security/main amd64 libicu67 amd64 67.1-7+deb11u1 [8624 kB]
Get:37 http://deb.debian.org/debian-security bullseye-security/main amd64 libxml2 amd64 2.9.10+dfsg-6.7+deb11u7 [693 kB]
Get:38 http://deb.debian.org/debian bullseye/main amd64 shared-mime-info amd64 2.0-1 [701 kB]
Get:39 http://deb.debian.org/debian bullseye/main amd64 xdg-user-dirs amd64 0.17-2 [53.8 kB]
Fetched 60.6 MB in 15s (4023 kB/s)
debconf: delaying package configuration, since apt-utils is not installed
Selecting previously unselected package libapparmor1:amd64.
(Reading database ... 6673 files and directories currently installed.)
Preparing to unpack .../00-libapparmor1_2.13.6-10_amd64.deb ...
Unpacking libapparmor1:amd64 (2.13.6-10) ...
Selecting previously unselected package libdbus-1-3:amd64.
Preparing to unpack .../01-libdbus-1-3_1.12.28-0+deb11u1_amd64.deb ...
Unpacking libdbus-1-3:amd64 (1.12.28-0+deb11u1) ...
Selecting previously unselected package libexpat1:amd64.
Preparing to unpack .../02-libexpat1_2.2.10-2+deb11u7_amd64.deb ...
Unpacking libexpat1:amd64 (2.2.10-2+deb11u7) ...
Selecting previously unselected package dbus.
Preparing to unpack .../03-dbus_1.12.28-0+deb11u1_amd64.deb ...
Unpacking dbus (1.12.28-0+deb11u1) ...
Selecting previously unselected package sensible-utils.
Preparing to unpack .../04-sensible-utils_0.0.14_all.deb ...
Unpacking sensible-utils (0.0.14) ...
Selecting previously unselected package openssl.
Preparing to unpack .../05-openssl_1.1.1w-0+deb11u3_amd64.deb ...
Unpacking openssl (1.1.1w-0+deb11u3) ...
Selecting previously unselected package ca-certificates.
Preparing to unpack .../06-ca-certificates_20210119_all.deb ...
Unpacking ca-certificates (20210119) ...
Selecting previously unselected package ucf.
Preparing to unpack .../07-ucf_3.0043+deb11u2_all.deb ...
Moving old data out of the way
Unpacking ucf (3.0043+deb11u2) ...
Selecting previously unselected package alsa-topology-conf.
Preparing to unpack .../08-alsa-topology-conf_1.2.4-1_all.deb ...
Unpacking alsa-topology-conf (1.2.4-1) ...
Selecting previously unselected package libasound2-data.
Preparing to unpack .../09-libasound2-data_1.2.4-1.1_all.deb ...
Unpacking libasound2-data (1.2.4-1.1) ...
Selecting previously unselected package libasound2:amd64.
Preparing to unpack .../10-libasound2_1.2.4-1.1_amd64.deb ...
Unpacking libasound2:amd64 (1.2.4-1.1) ...
Selecting previously unselected package alsa-ucm-conf.
Preparing to unpack .../11-alsa-ucm-conf_1.2.4-2_all.deb ...
Unpacking alsa-ucm-conf (1.2.4-2) ...
Selecting previously unselected package java-common.
Preparing to unpack .../12-java-common_0.72_all.deb ...
Unpacking java-common (0.72) ...
Selecting previously unselected package libavahi-common-data:amd64.
Preparing to unpack .../13-libavahi-common-data_0.8-5+deb11u3_amd64.deb ...
Unpacking libavahi-common-data:amd64 (0.8-5+deb11u3) ...
Selecting previously unselected package libavahi-common3:amd64.
Preparing to unpack .../14-libavahi-common3_0.8-5+deb11u3_amd64.deb ...
Unpacking libavahi-common3:amd64 (0.8-5+deb11u3) ...
Selecting previously unselected package libavahi-client3:amd64.
Preparing to unpack .../15-libavahi-client3_0.8-5+deb11u3_amd64.deb ...
Unpacking libavahi-client3:amd64 (0.8-5+deb11u3) ...
Selecting previously unselected package libcups2:amd64.
Preparing to unpack .../16-libcups2_2.3.3op2-3+deb11u9_amd64.deb ...
Unpacking libcups2:amd64 (2.3.3op2-3+deb11u9) ...
Selecting previously unselected package libbrotli1:amd64.
Preparing to unpack .../17-libbrotli1_1.0.9-2+b2_amd64.deb ...
Unpacking libbrotli1:amd64 (1.0.9-2+b2) ...
Selecting previously unselected package libpng16-16:amd64.
Preparing to unpack .../18-libpng16-16_1.6.37-3_amd64.deb ...
Unpacking libpng16-16:amd64 (1.6.37-3) ...
Selecting previously unselected package libfreetype6:amd64.
Preparing to unpack .../19-libfreetype6_2.10.4+dfsg-1+deb11u2_amd64.deb ...
Unpacking libfreetype6:amd64 (2.10.4+dfsg-1+deb11u2) ...
Selecting previously unselected package fonts-dejavu-core.
Preparing to unpack .../20-fonts-dejavu-core_2.37-2_all.deb ...
Unpacking fonts-dejavu-core (2.37-2) ...
Selecting previously unselected package fontconfig-config.
Preparing to unpack .../21-fontconfig-config_2.13.1-4.2_all.deb ...
Unpacking fontconfig-config (2.13.1-4.2) ...
Selecting previously unselected package libfontconfig1:amd64.
Preparing to unpack .../22-libfontconfig1_2.13.1-4.2_amd64.deb ...
Unpacking libfontconfig1:amd64 (2.13.1-4.2) ...
Selecting previously unselected package liblcms2-2:amd64.
Preparing to unpack .../23-liblcms2-2_2.12~rc1-2_amd64.deb ...
Unpacking liblcms2-2:amd64 (2.12~rc1-2) ...
Selecting previously unselected package libjpeg62-turbo:amd64.
Preparing to unpack .../24-libjpeg62-turbo_1%3a2.0.6-4_amd64.deb ...
Unpacking libjpeg62-turbo:amd64 (1:2.0.6-4) ...
Selecting previously unselected package libnspr4:amd64.
Preparing to unpack .../25-libnspr4_2%3a4.29-1_amd64.deb ...
Unpacking libnspr4:amd64 (2:4.29-1) ...
Selecting previously unselected package libsqlite3-0:amd64.
Preparing to unpack .../26-libsqlite3-0_3.34.1-3+deb11u1_amd64.deb ...
Unpacking libsqlite3-0:amd64 (3.34.1-3+deb11u1) ...
Selecting previously unselected package libnss3:amd64.
Preparing to unpack .../27-libnss3_2%3a3.61-1+deb11u4_amd64.deb ...
Unpacking libnss3:amd64 (2:3.61-1+deb11u4) ...
Selecting previously unselected package libglib2.0-0:amd64.
Preparing to unpack .../28-libglib2.0-0_2.66.8-1+deb11u6_amd64.deb ...
Unpacking libglib2.0-0:amd64 (2.66.8-1+deb11u6) ...
Selecting previously unselected package libgraphite2-3:amd64.
Preparing to unpack .../29-libgraphite2-3_1.3.14-1_amd64.deb ...
Unpacking libgraphite2-3:amd64 (1.3.14-1) ...
Selecting previously unselected package libharfbuzz0b:amd64.
Preparing to unpack .../30-libharfbuzz0b_2.7.4-1_amd64.deb ...
Unpacking libharfbuzz0b:amd64 (2.7.4-1) ...
Selecting previously unselected package libpcsclite1:amd64.
Preparing to unpack .../31-libpcsclite1_1.9.1-1_amd64.deb ...
Unpacking libpcsclite1:amd64 (1.9.1-1) ...
Selecting previously unselected package openjdk-11-jre-headless:amd64.
Preparing to unpack .../32-openjdk-11-jre-headless_11.0.27+6-1~deb11u1_amd64.deb ...
Unpacking openjdk-11-jre-headless:amd64 (11.0.27+6-1~deb11u1) ...
Selecting previously unselected package ca-certificates-java.
Preparing to unpack .../33-ca-certificates-java_20190909+deb11u1_all.deb ...
Unpacking ca-certificates-java (20190909+deb11u1) ...
Selecting previously unselected package libglib2.0-data.
Preparing to unpack .../34-libglib2.0-data_2.66.8-1+deb11u6_all.deb ...
Unpacking libglib2.0-data (2.66.8-1+deb11u6) ...
Selecting previously unselected package libicu67:amd64.
Preparing to unpack .../35-libicu67_67.1-7+deb11u1_amd64.deb ...
Unpacking libicu67:amd64 (67.1-7+deb11u1) ...
Selecting previously unselected package libxml2:amd64.
Preparing to unpack .../36-libxml2_2.9.10+dfsg-6.7+deb11u7_amd64.deb ...
Unpacking libxml2:amd64 (2.9.10+dfsg-6.7+deb11u7) ...
Selecting previously unselected package shared-mime-info.
Preparing to unpack .../37-shared-mime-info_2.0-1_amd64.deb ...
Unpacking shared-mime-info (2.0-1) ...
Selecting previously unselected package xdg-user-dirs.
Preparing to unpack .../38-xdg-user-dirs_0.17-2_amd64.deb ...
Unpacking xdg-user-dirs (0.17-2) ...
Setting up libexpat1:amd64 (2.2.10-2+deb11u7) ...
Setting up libgraphite2-3:amd64 (1.3.14-1) ...
Setting up liblcms2-2:amd64 (2.12~rc1-2) ...
Setting up libapparmor1:amd64 (2.13.6-10) ...
Setting up java-common (0.72) ...
Setting up libicu67:amd64 (67.1-7+deb11u1) ...
Setting up xdg-user-dirs (0.17-2) ...
Setting up libglib2.0-0:amd64 (2.66.8-1+deb11u6) ...
No schema files found: doing nothing.
Setting up libbrotli1:amd64 (1.0.9-2+b2) ...
Setting up libsqlite3-0:amd64 (3.34.1-3+deb11u1) ...
Setting up libasound2-data (1.2.4-1.1) ...
Setting up libglib2.0-data (2.66.8-1+deb11u6) ...
Setting up libjpeg62-turbo:amd64 (1:2.0.6-4) ...
Setting up libnspr4:amd64 (2:4.29-1) ...
Setting up libavahi-common-data:amd64 (0.8-5+deb11u3) ...
Setting up libdbus-1-3:amd64 (1.12.28-0+deb11u1) ...
Setting up dbus (1.12.28-0+deb11u1) ...
invoke-rc.d: could not determine current runlevel
invoke-rc.d: policy-rc.d denied execution of start.
Setting up libpng16-16:amd64 (1.6.37-3) ...
Setting up fonts-dejavu-core (2.37-2) ...
Setting up libpcsclite1:amd64 (1.9.1-1) ...
Setting up alsa-topology-conf (1.2.4-1) ...
Setting up sensible-utils (0.0.14) ...
Setting up libasound2:amd64 (1.2.4-1.1) ...
Setting up openssl (1.1.1w-0+deb11u3) ...
Setting up libxml2:amd64 (2.9.10+dfsg-6.7+deb11u7) ...
Setting up alsa-ucm-conf (1.2.4-2) ...
Setting up libavahi-common3:amd64 (0.8-5+deb11u3) ...
Setting up libnss3:amd64 (2:3.61-1+deb11u4) ...
Setting up ca-certificates (20210119) ...
debconf: unable to initialize frontend: Dialog
debconf: (No usable dialog-like program is installed, so the dialog based frontend cannot be used. at /usr/share/perl5/Debconf/FrontEnd/Dialog.pm line 78.)
debconf: falling back to frontend: Readline
debconf: unable to initialize frontend: Readline
debconf: (Can't locate Term/ReadLine.pm in @INC (you may need to install the Term::ReadLine module) (@INC contains: /etc/perl /usr/local/lib/x86_64-linux-gnu/perl/5.32.1 /usr/local/share/perl/5.32.1 /usr/lib/x86_64-linux-gnu/perl5/5.32 /usr/share/perl5 /usr/lib/x86_64-linux-gnu/perl-base /usr/lib/x86_64-linux-gnu/perl/5.32 /usr/share/perl/5.32 /usr/local/lib/site_perl) at /usr/share/perl5/Debconf/FrontEnd/Readline.pm line 7.)
debconf: falling back to frontend: Teletype
Updating certificates in /etc/ssl/certs...
129 added, 0 removed; done.
Setting up libfreetype6:amd64 (2.10.4+dfsg-1+deb11u2) ...
Setting up shared-mime-info (2.0-1) ...
Setting up ucf (3.0043+deb11u2) ...
debconf: unable to initialize frontend: Dialog
debconf: (No usable dialog-like program is installed, so the dialog based frontend cannot be used. at /usr/share/perl5/Debconf/FrontEnd/Dialog.pm line 78.)
debconf: falling back to frontend: Readline
debconf: unable to initialize frontend: Readline
debconf: (Can't locate Term/ReadLine.pm in @INC (you may need to install the Term::ReadLine module) (@INC contains: /etc/perl /usr/local/lib/x86_64-linux-gnu/perl/5.32.1 /usr/local/share/perl/5.32.1 /usr/lib/x86_64-linux-gnu/perl5/5.32 /usr/share/perl5 /usr/lib/x86_64-linux-gnu/perl-base /usr/lib/x86_64-linux-gnu/perl/5.32 /usr/share/perl/5.32 /usr/local/lib/site_perl) at /usr/share/perl5/Debconf/FrontEnd/Readline.pm line 7.)
debconf: falling back to frontend: Teletype
Setting up libharfbuzz0b:amd64 (2.7.4-1) ...
Setting up libavahi-client3:amd64 (0.8-5+deb11u3) ...
Setting up fontconfig-config (2.13.1-4.2) ...
debconf: unable to initialize frontend: Dialog
debconf: (No usable dialog-like program is installed, so the dialog based frontend cannot be used. at /usr/share/perl5/Debconf/FrontEnd/Dialog.pm line 78.)
debconf: falling back to frontend: Readline
debconf: unable to initialize frontend: Readline
debconf: (Can't locate Term/ReadLine.pm in @INC (you may need to install the Term::ReadLine module) (@INC contains: /etc/perl /usr/local/lib/x86_64-linux-gnu/perl/5.32.1 /usr/local/share/perl/5.32.1 /usr/lib/x86_64-linux-gnu/perl5/5.32 /usr/share/perl5 /usr/lib/x86_64-linux-gnu/perl-base /usr/lib/x86_64-linux-gnu/perl/5.32 /usr/share/perl/5.32 /usr/local/lib/site_perl) at /usr/share/perl5/Debconf/FrontEnd/Readline.pm line 7.)
debconf: falling back to frontend: Teletype
Setting up libcups2:amd64 (2.3.3op2-3+deb11u9) ...
Setting up libfontconfig1:amd64 (2.13.1-4.2) ...
Setting up ca-certificates-java (20190909+deb11u1) ...
'/etc/java-11-openjdk/security/java.security.dpkg-new' -> '/etc/java-11-openjdk/security/java.security'
head: cannot open '/etc/ssl/certs/java/cacerts' for reading: No such file or directory
Adding debian:Hongkong_Post_Root_CA_1.pem
Adding debian:TeliaSonera_Root_CA_v1.pem
Adding debian:D-TRUST_Root_Class_3_CA_2_EV_2009.pem
Adding debian:GlobalSign_Root_CA_-_R2.pem
Adding debian:COMODO_ECC_Certification_Authority.pem
Adding debian:VeriSign_Universal_Root_Certification_Authority.pem
Adding debian:Certum_Trusted_Network_CA_2.pem
Adding debian:Comodo_AAA_Services_root.pem
Adding debian:Hellenic_Academic_and_Research_Institutions_RootCA_2011.pem
Adding debian:USERTrust_ECC_Certification_Authority.pem
Adding debian:AffirmTrust_Networking.pem
Adding debian:Trustwave_Global_ECC_P256_Certification_Authority.pem
Adding debian:TWCA_Global_Root_CA.pem
Adding debian:TWCA_Root_Certification_Authority.pem
Adding debian:SecureTrust_CA.pem
Adding debian:ISRG_Root_X1.pem
Adding debian:COMODO_RSA_Certification_Authority.pem
Adding debian:QuoVadis_Root_CA.pem
Adding debian:Certigna.pem
Adding debian:DigiCert_Global_Root_CA.pem
Adding debian:Trustwave_Global_ECC_P384_Certification_Authority.pem
Adding debian:EC-ACC.pem
Adding debian:COMODO_Certification_Authority.pem
Adding debian:ACCVRAIZ1.pem
Adding debian:certSIGN_Root_CA_G2.pem
Adding debian:DigiCert_High_Assurance_EV_Root_CA.pem
Adding debian:QuoVadis_Root_CA_2.pem
Adding debian:Hongkong_Post_Root_CA_3.pem
Adding debian:GTS_Root_R4.pem
Adding debian:Starfield_Services_Root_Certificate_Authority_-_G2.pem
Adding debian:Amazon_Root_CA_3.pem
Adding debian:Baltimore_CyberTrust_Root.pem
Adding debian:Entrust_Root_Certification_Authority_-_G4.pem
Adding debian:Go_Daddy_Class_2_CA.pem
Adding debian:Amazon_Root_CA_4.pem
Adding debian:emSign_ECC_Root_CA_-_G3.pem
Adding debian:DigiCert_Global_Root_G3.pem
Adding debian:Starfield_Class_2_CA.pem
Adding debian:Chambers_of_Commerce_Root_-_2008.pem
Adding debian:Actalis_Authentication_Root_CA.pem
Adding debian:GTS_Root_R2.pem
Adding debian:Atos_TrustedRoot_2011.pem
Adding debian:NetLock_Arany_=Class_Gold=_Főtanúsítvány.pem
Adding debian:SSL.com_EV_Root_Certification_Authority_RSA_R2.pem
Adding debian:DigiCert_Assured_ID_Root_G2.pem
Adding debian:Staat_der_Nederlanden_Root_CA_-_G3.pem
Adding debian:certSIGN_ROOT_CA.pem
Adding debian:SZAFIR_ROOT_CA2.pem
Adding debian:SecureSign_RootCA11.pem
Adding debian:Izenpe.com.pem
Adding debian:emSign_ECC_Root_CA_-_C3.pem
Adding debian:Buypass_Class_3_Root_CA.pem
Adding debian:AffirmTrust_Commercial.pem
Adding debian:DigiCert_Global_Root_G2.pem
Adding debian:TUBITAK_Kamu_SM_SSL_Kok_Sertifikasi_-_Surum_1.pem
Adding debian:Secure_Global_CA.pem
Adding debian:QuoVadis_Root_CA_2_G3.pem
Adding debian:SwissSign_Silver_CA_-_G2.pem
Adding debian:SSL.com_EV_Root_Certification_Authority_ECC.pem
Adding debian:Trustwave_Global_Certification_Authority.pem
Adding debian:UCA_Global_G2_Root.pem
Adding debian:TrustCor_ECA-1.pem
Adding debian:TrustCor_RootCert_CA-1.pem
Adding debian:emSign_Root_CA_-_G1.pem
Adding debian:Go_Daddy_Root_Certificate_Authority_-_G2.pem
Adding debian:Certigna_Root_CA.pem
Adding debian:GlobalSign_Root_CA_-_R6.pem
Adding debian:Entrust.net_Premium_2048_Secure_Server_CA.pem
Adding debian:UCA_Extended_Validation_Root.pem
Adding debian:Entrust_Root_Certification_Authority_-_EC1.pem
Adding debian:Autoridad_de_Certificacion_Firmaprofesional_CIF_A62634068.pem
Adding debian:ePKI_Root_Certification_Authority.pem
Adding debian:e-Szigno_Root_CA_2017.pem
Adding debian:Microsoft_RSA_Root_Certificate_Authority_2017.pem
Adding debian:Buypass_Class_2_Root_CA.pem
Adding debian:CFCA_EV_ROOT.pem
Adding debian:NAVER_Global_Root_Certification_Authority.pem
Adding debian:GlobalSign_ECC_Root_CA_-_R4.pem
Adding debian:DigiCert_Trusted_Root_G4.pem
Adding debian:GTS_Root_R3.pem
Adding debian:T-TeleSec_GlobalRoot_Class_3.pem
Adding debian:Amazon_Root_CA_1.pem
Adding debian:AffirmTrust_Premium_ECC.pem
Adding debian:E-Tugra_Certification_Authority.pem
Adding debian:GlobalSign_Root_CA.pem
Adding debian:Security_Communication_Root_CA.pem
Adding debian:USERTrust_RSA_Certification_Authority.pem
Adding debian:QuoVadis_Root_CA_3.pem
Adding debian:Global_Chambersign_Root_-_2008.pem
Adding debian:Entrust_Root_Certification_Authority.pem
Adding debian:emSign_Root_CA_-_C1.pem
Adding debian:Trustis_FPS_Root_CA.pem
Adding debian:QuoVadis_Root_CA_3_G3.pem
Adding debian:AC_RAIZ_FNMT-RCM.pem
Adding debian:SSL.com_Root_Certification_Authority_RSA.pem
Adding debian:DST_Root_CA_X3.pem
Adding debian:IdenTrust_Commercial_Root_CA_1.pem
Adding debian:GlobalSign_ECC_Root_CA_-_R5.pem
Adding debian:GDCA_TrustAUTH_R5_ROOT.pem
Adding debian:SwissSign_Gold_CA_-_G2.pem
Adding debian:Entrust_Root_Certification_Authority_-_G2.pem
Adding debian:Certum_Trusted_Network_CA.pem
Adding debian:SSL.com_Root_Certification_Authority_ECC.pem
Adding debian:AffirmTrust_Premium.pem
Adding debian:D-TRUST_Root_Class_3_CA_2_2009.pem
Adding debian:CA_Disig_Root_R2.pem
Adding debian:Hellenic_Academic_and_Research_Institutions_RootCA_2015.pem
Adding debian:QuoVadis_Root_CA_1_G3.pem
Adding debian:Staat_der_Nederlanden_EV_Root_CA.pem
Adding debian:Network_Solutions_Certificate_Authority.pem
Adding debian:GeoTrust_Primary_Certification_Authority_-_G2.pem
Adding debian:TrustCor_RootCert_CA-2.pem
Adding debian:GlobalSign_Root_CA_-_R3.pem
Adding debian:XRamp_Global_CA_Root.pem
Adding debian:T-TeleSec_GlobalRoot_Class_2.pem
Adding debian:Microsec_e-Szigno_Root_CA_2009.pem
Adding debian:Amazon_Root_CA_2.pem
Adding debian:OISTE_WISeKey_Global_Root_GC_CA.pem
Adding debian:IdenTrust_Public_Sector_Root_CA_1.pem
Adding debian:Starfield_Root_Certificate_Authority_-_G2.pem
Adding debian:Sonera_Class_2_Root_CA.pem
Adding debian:Microsoft_ECC_Root_Certificate_Authority_2017.pem
Adding debian:DigiCert_Assured_ID_Root_CA.pem
Adding debian:Cybertrust_Global_Root.pem
Adding debian:GTS_Root_R1.pem
Adding debian:OISTE_WISeKey_Global_Root_GB_CA.pem
Adding debian:Hellenic_Academic_and_Research_Institutions_ECC_RootCA_2015.pem
Adding debian:Security_Communication_RootCA2.pem
Adding debian:DigiCert_Assured_ID_Root_G3.pem
done.
Processing triggers for libc-bin (2.31-13+deb11u13) ...
Processing triggers for ca-certificates (20210119) ...
Updating certificates in /etc/ssl/certs...
0 added, 0 removed; done.
Running hooks in /etc/ca-certificates/update.d...

done.
done.
Setting up openjdk-11-jre-headless:amd64 (11.0.27+6-1~deb11u1) ...
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/java to provide /usr/bin/java (java) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/jjs to provide /usr/bin/jjs (jjs) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/keytool to provide /usr/bin/keytool (keytool) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/rmid to provide /usr/bin/rmid (rmid) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/rmiregistry to provide /usr/bin/rmiregistry (rmiregistry) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/pack200 to provide /usr/bin/pack200 (pack200) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/bin/unpack200 to provide /usr/bin/unpack200 (unpack200) in auto mode
update-alternatives: using /usr/lib/jvm/java-11-openjdk-amd64/lib/jexec to provide /usr/bin/jexec (jexec) in auto mode
root@79cf31c19188:/# cd mnt/cplex_stuff/
root@79cf31c19188:/mnt/cplex_stuff# ls
cplex_studio2212.linux_x86_64.bin
root@79cf31c19188:/mnt/cplex_stuff# ./cplex_studio2212.linux_x86_64.bin
Preparing to install
Extracting the JRE from the installer archive...
Unpacking the JRE...
Extracting the installation resources from the installer archive...
Configuring the installer for this system's environment...

Launching installer...

===============================================================================
Choose Locale...
----------------

    1- Deutsch
  ->2- English
    3- Espa?ol
    4- Fran?ais
    5- Italiano
    6- Portugu?s  (Brasil)

CHOOSE LOCALE BY NUMBER: 2
===============================================================================
IBM ILOG CPLEX Optimization Studio 22.1.2        (created with InstallAnywhere)
-------------------------------------------------------------------------------

Preparing CONSOLE Mode Installation...




===============================================================================
Introduction
------------

InstallAnywhere will guide you through the installation of IBM ILOG CPLEX
Optimization Studio 22.1.2.

It is strongly recommended that you quit all programs before continuing with
this installation.

Respond to each prompt to proceed to the next step in the installation.  If
you want to change something on a previous step, type 'back'.

You may cancel this installation at any time by typing 'quit'.

PRESS <ENTER> TO CONTINUE:



===============================================================================




    LICENSE INFORMATION

    The Programs listed below are licensed under the following License
    Information terms and conditions in addition to the Program license
    terms previously agreed to by Client and IBM. If Client does not have
    previously agreed to license terms in effect for the Program, the
    International Program License Agreement (i125-3301-15) applies.

    Program Name (Program Number):
    IBM ILOG OPL CPLEX Developer Edition 22.1.2 (5724-Y54)
    IBM ILOG CPLEX MILP add-on Desktop Edition 22.1.2 (5724-Y48)
    IBM ILOG CPLEX Desktop Edition 22.1.2 (5724-Y48)
    IBM ILOG CPLEX LP Desktop Edition 22.1.2 (5724-Y48)
    IBM ILOG CP Optimizer Desktop Edition 22.1.2 (5724-Y49)
    IBM ILOG CPLEX For Non Production 22.1.2 (5724-Y48)
    IBM ILOG CPLEX MILP add-on Server Edition 22.1.2 (5724-Y48)
    IBM ILOG CPLEX Server Edition 22.1.2 (5724-Y48)
    IBM ILOG CPLEX LP Server Edition 22.1.2 (5724-Y48)
    IBM ILOG CP Optimizer Server Edition 22.1.2 (5724-Y49)

Press Enter to continue viewing the license agreement, or enter "1" to
   accept the agreement, "2" to decline it, "3" to print it, or "99" to go back
   to the previous screen.: 1




===============================================================================
Choose Install Folder
---------------------

Where would you like to install?

  Default Install Folder: /opt/ibm/ILOG/CPLEX_Studio2212

ENTER AN ABSOLUTE PATH, OR PRESS <ENTER> TO ACCEPT THE DEFAULT
      :



===============================================================================
Ready To Install
----------------

InstallAnywhere is now ready to install IBM ILOG CPLEX Optimization Studio
22.1.2 onto your system at the following location:

   /opt/ibm/ILOG/CPLEX_Studio2212

PRESS <ENTER> TO INSTALL:



===============================================================================
Pre-Installation Summary
------------------------

Please Review the Following Before Continuing:

Product Name:
    IBM ILOG CPLEX Optimization Studio 22.1.2

Install Folder:
    /opt/ibm/ILOG/CPLEX_Studio2212

Product Version
    22.1.2

Disk Space Information (for Installation Target):
    Required:    1,225,336,428 Bytes
    Available: 983,726,784,512 Bytes

PRESS <ENTER> TO CONTINUE:



===============================================================================
Installing...
-------------

 [==================|==================|==================|==================]
 [------------------|------------------|------------------|-----------------


===============================================================================
Prepare env

Please Wait
-----------
-]



===============================================================================


Please Wait
-----------



===============================================================================


Please Wait
-----------



===============================================================================
Installation Complete
---------------------

IBM ILOG CPLEX Optimization Studio 22.1.2 has been successfully installed to:

   /opt/ibm/ILOG/CPLEX_Studio2212

PRESS <ENTER> TO EXIT THE INSTALLER:
root@79cf31c19188:/mnt/cplex_stuff# cd ..
root@79cf31c19188:/mnt# cd ..
root@79cf31c19188:/# cd opt/
root@79cf31c19188:/opt# ls
ibm
root@79cf31c19188:/opt# tar czf /mnt/cplex.tar.gz -C /opt ibm
root@79cf31c19188:/opt# exit
exit
PS C:\Users\agarg\research\Microbiome-GEM-Modeling>
