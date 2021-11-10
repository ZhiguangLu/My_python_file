from email.mime.text import MIMEText
from email.header import Header
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
import zipfile
import smtplib
import os
import datetime


class Email:
    """
    from_mail: 发送者的邮箱 (字符串形式)
    mail_pass： 发送者邮箱的授权码 (字符串形式)注意：授权码不是登录密码，需要自己在邮箱开启，可以百度邮箱授权码
    to_mail： 收件人的邮箱 (列表类似，[字符串1,字符串2,...],每个字符串都是一个邮箱)
    dir_path： 需要压缩的文件夹路径 (字符串形式，路径是用这种/斜杠)
    outFullName： 压缩文件保存的路径+压缩文件名称(如Name.rar) (字符串形式)
    content: 邮件正文内容 (字符串形式)
    topic: 邮件主题 (字符串形式)
    """

    def __init__(self, from_mail, mail_pass, to_mail, dir_path, outFullName, content, topic):
        self.from_mail = from_mail
        self.mail_pass = mail_pass
        self.to_mail = to_mail
        self.dir_path = dir_path
        self.outFullName = outFullName
        self.content = content
        self.topic = topic
        self.rar_dir()
        self.send_email()

    def rar_dir(self):
        """
        压缩指定文件夹
        """
        testcase_zip = zipfile.ZipFile(self.outFullName, 'w', zipfile.ZIP_DEFLATED)
        for path, dir_names, file_names in os.walk(self.dir_path):
            for filename in file_names:
                testcase_zip.write(os.path.join(path, filename))
        testcase_zip.close()
        print("文件夹压缩成功")

    def set_content(self):
        path = self.outFullName
        name_of_rar = path[path.rfind('/') + 1:]
        if 'qq' in self.from_mail:
            smtp_server = 'smtp.qq.com'
        else:
            smtp_server = 'smtp.163.com'
        msg = MIMEMultipart()
        msg['From'] = self.from_mail
        msg['To'] = ','.join(self.to_mail)
        msg['Subject'] = Header(self.topic, 'utf-8').encode()
        msg.attach(MIMEText(self.content, 'plain', 'utf-8'))

        with open(name_of_rar, 'rb') as f:
            mime = MIMEBase('rar', 'rar', filename=name_of_rar)
            mime.add_header('Content-Disposition', 'attachment', filename=name_of_rar)
            mime.set_payload(f.read())
        encoders.encode_base64(mime)
        msg.attach(mime)
        return msg, smtp_server

    def send_email(self):
        msg, smep_server = self.set_content()
        try:
            s = smtplib.SMTP()
            s.connect(smep_server, "25")
            s.login(self.from_mail, self.mail_pass)
            s.sendmail(self.from_mail, self.to_mail, msg.as_string())
            s.quit()
            print('邮件发送成功！')
        except smtplib.SMTPException as e:
            print("Error: %s" % e)


f_e = '13246457757@163.com'
e_p = 'IBUGSNANUVEDGZKD'
t_e = ['497025048@qq.com']
d_p = 'C:/Users/Young/Desktop/picture'
o_p = 'E:/Pycharm/Python_text/exercise/test.rar'
c_t = '时间{}程序已经跑完，现在将所得数据以邮件形式发送给你，请接收！！！'.format(str(datetime.datetime.today()))
t_p = 'fisher_information_test'

Email(f_e, e_p, t_e, d_p, o_p, c_t, t_p)
