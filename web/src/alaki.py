from threading import Thread
from time import sleep # for time.sleep


class MYCLASS:
    def mymain(self):
        thread = Thread(target=self.func);
        thread.start()
        for i in range(5):
            print 'mymain_', i
            sleep(1)
        
        thread.join()
    
    def func(self):
        for i in range(5):
            print 'func_', i
            sleep(1)

if __name__ == '__main__':
    c = MYCLASS()
    c.mymain()