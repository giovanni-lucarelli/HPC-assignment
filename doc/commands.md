# Comandi utili

Copiare un file dal proprio computer a un server remoto usando `scp`:

```bash

scp -r /home/giovanni/Desktop/HPC-orfeo gioluc@orfeo:/u/ipauser/gioluc

```

Copiare un file da un server remoto al proprio computer usando `scp`:

```bash

scp -r gioluc@orfeo:/u/ipauser/gioluc/HPC-orfeo /home/giovanni/Desktop/HPC-orfeo

```

> Nota: entrambi i comandi si eseguono da terminale sul proprio computer (non sul server remoto!).
