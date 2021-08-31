# ローカル環境で作成したHPをremote repositoryへ反映させる(pushする)方法
Git Bashでローカル環境の「v2mv02e2ax.github.io」のディレクトリへ移動して以下のコードを打ち込む。
```bash
git add .
git commit -m "first commit"
git push origin main
```

# remote repositoryで変更したHPをローカル環境へ反映させる（pullする）方法
remote repositoryでファイルを変更してpullをせずに、ローカル環境にpushするとconflictするので注意。
Git BashでLocal環境の「v2mv02e2ax.github.io」のディレクトリへ移動して以下のコードを打ち込む。
```bash
git pull
```
